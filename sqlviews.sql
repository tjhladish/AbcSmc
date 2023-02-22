
/* TODO: should jobs also be some sort of view? */
CREATE TABLE IF NOT EXISTS job (
    serial INTEGER PRIMARY KEY ASC, -- all runs must have a serial
    smcSet INTEGER,                 -- smcSet only required when doing fitting
    particleIdx INTEGER,            -- particleIdx only required when doing scenario analysis
    startTime INTEGER,              -- start time of run; initially NULL, then set at time of run
    duration INTEGER,               -- duration of run; initially NULL, then set at time of run
    status TEXT DEFAULT 'Q',        -- status of run; initially 'Q' for queued, then 'D' for done
    posterior INTEGER DEFAULT -1,   -- posterior rank; defaults to -1 (not in posterior), potentially updated by processing
    attempts INTEGER DEFAULT 0      -- number of times this job has been attempted; initially 0, then incremented on "checkout"
);

/*
CREATE TABLE IF NOT EXISTS job_data (
    serial INTEGER PRIMARY KEY ASC,
    startTime INTEGER, duration INTEGER,
    attempts INTEGER DEFAULT 0,
    priority INTEGER
);

CREATE TABLE IF NOT EXISTS set_data (
    serial INTEGER PRIMARY KEY ASC,
    setId INTEGER DEFAULT -1,
    rank INTEGER DEFAULT -1
);

CREATE VIEW IF NOT EXISTS job AS
    SELECT serial, COALESCE(setId, -1) AS smcSet, 
    startTime, duration,
    status,
    rank AS posterior, attempts
    FROM job_data
    LEFT JOIN set_data USING (serial)
    LEFT JOIN met_vals USING (serial)
    WHERE metIdx = 1;
*/

CREATE INDEX IF NOT EXISTS idx1 ON job (status, attempts);

/* random seed container; "seed" is a special kind of parameter */
CREATE TABLE IF NOT EXISTS seeds (
    serial INTEGER PRIMARY KEY ASC,
    seed BLOB
);

/* normal-formed parameters table. links run serial + parameter id => value of that parameter */
CREATE TABLE IF NOT EXISTS par_vals (
    serial INTEGER NOT NULL,
    parIdx INTEGER NOT NULL,
    value REAL,
    PRIMARY KEY (serial, parIdx)
);

/* normal-formed untransformed parameters table. links run serial + parameter id => untrans value of that parameter */
/* always created, but can be left empty and everything else will just work. */
CREATE TABLE IF NOT EXISTS upar_vals (
    serial INTEGER NOT NULL,
    parIdx INTEGER NOT NULL,
    uvalue REAL,
    PRIMARY KEY (serial, parIdx)
);

/* meta-data on parameters; must be dynamically filled when setting up db */
/* TODO: should this include other meta information, e.g. if a parameter is untransformed? */
CREATE TABLE IF NOT EXISTS par_name (
    parIdx INTEGER PRIMARY KEY ASC,
    name TEXT NOT NULL,
    long_name TEXT
);

/* normal-formed metrics table. links run serial + metric id => value of that metric */
CREATE TABLE IF NOT EXISTS met_vals (
    serial INTEGER NOT NULL,
    metIdx INTEGER NOT NULL,
    value REAL,
    PRIMARY KEY (serial, metIdx)
);

/* meta-data on metrics; must be dynamically filled when setting up db */
/* TODO: should this include other meta information, e.g. the observed value? */
CREATE TABLE IF NOT EXISTS met_name (
    metIdx INTEGER PRIMARY KEY ASC,
    name TEXT NOT NULL,
    long_name TEXT
);

/* SQLite supports VIEWs - these are essentially predefined queries, that are dynamically
*  evaluated, and function as-if they were TABLEs.
*
*  For `AbcSmc` SQLite-based storage, we define several views which encapsulate queries that
*  represent common operations as part of the analysis flow.
*/

/* ================== SMC Set Level Operations & Posterior Queries ===================== */

/* add the "smc set summary" view: this provides a convenient summary of the stages of analysis
*  This view is also used as the starting point for other views dedicated to gathering smc sets
*  by analysis stage.
*/
/* smc sets where total == run AND postsize > 0 are considered complete */
/* smc sets where total == run AND postsize == 0 are considered TODO */
/* smc sets where total > run => still has work to be simulated before being ready to process */
/* smcSet represents ascending order of runs, so the "last" set meeting some criteria is MAX(smcSet)
*  and the "next" set meeting some criteria is MIN(smcSet)
*/
/* TODO: if configuration data were available in db, could assess if the
*  "right" number of jobs were run, and if the "right" number of jobs have
*  posterior > -1. Right now have to rely on the algorithmic guarantee that
*  writing posterior ranks matches the configuration file + writes them all
*  as a transaction.
*/
CREATE VIEW IF NOT EXISTS smc_summary AS
    SELECT smcSet, MIN(serial) AS start_serial, MAX(serial) AS end_serial,
        COUNT(*) AS total, SUM(CASE WHEN status == 'D' THEN 1 ELSE 0 END) AS run,
        SUM(CASE WHEN posterior != -1 THEN 1 ELSE 0 END) AS postsize
    FROM job
    GROUP BY smcSet ORDER BY smcSet;

/* add the "complete posterior" smc waves view: has all the information
*  needed to pull out parameters, metrics, by smcSet. Definitionally
*  ignores elements *not* in the posterior and makes no reference to job status
*  etc, because that's defined as done by this view.
*/
CREATE VIEW IF NOT EXISTS post_done AS
    SELECT serial, smcSet, posterior
    FROM job
    JOIN (SELECT smcSet FROM smc_summary WHERE total == run AND postsize != 0)
    USING (smcSet)
    ORDER BY smcSet, posterior;

/* add the "last complete posterior" view: finds the max value
*  and then filters on that (via join - cannot directly use WHERE smcSet == MAX(smcSet))
*/
CREATE VIEW IF NOT EXISTS post_last AS
    SELECT * FROM post_done JOIN (
        SELECT MAX(smcSet) AS smcSet FROM post_done
    ) USING (smcSet)
    ORDER BY posterior;

/* add the "still to process" view */
/* NB: not aware of uses were there are multiple done-but-not processed waves, but
*  can imagine a situation where that might be useful? Hence still introducing
*  the all (post_todo) vs one (post_next) distinction
*/
CREATE VIEW IF NOT EXISTS post_todo AS
    SELECT serial, smcSet
    FROM job
    JOIN (SELECT smcSet FROM smc_summary WHERE total == run AND posterior == 0)
    USING (smcSet)
    ORDER BY smcSet, posterior;

/* add the "next to process" view;
*  note parallel structure of post_todo/post_next vs post_done/post_last
*/
CREATE VIEW IF NOT EXISTS post_next AS
    SELECT * FROM post_todo JOIN (
        SELECT MIN(smcSet) AS smcSet FROM post_todo
    ) USING (smcSet)
    ORDER BY posterior;

/* add the "still to simulate" view
*  TODO: how to deal with prioritization of work?
*/
CREATE VIEW IF NOT EXISTS work_open AS
    SELECT smcSet, priority, serial
    FROM job
    JOIN (SELECT serial, CASE status WHEN 'Q' THEN attempts ELSE attempts + 100 END AS priority FROM job) USING (serial)
    WHERE status NOT IN ('D', 'P')
    ORDER BY smcSet, priority, serial;

/* add the "next to simulate" view */
CREATE VIEW IF NOT EXISTS work_next AS
    SELECT serial FROM work_open JOIN (
        SELECT MIN(smcSet) AS smcSet FROM work_open
    ) USING (smcSet)
    ORDER BY priority, serial;