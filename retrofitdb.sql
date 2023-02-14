/* SQLite supports VIEWs - these are essentially predefined queries, that are dynamically
*  evaluated, and function as-if they were TABLEs.
*
*  For `AbcSmc` SQLite-based storage, we define several views which encapsulate queries that
*  represent common operations as part of the analysis flow.
*/

/* ================== SMC Set Level Operations & Posterior Queries ===================== */

/* add the "smc set summary" view */
CREATE VIEW IF NOT EXISTS smc_summary AS
    SELECT smcSet, MIN(serial) AS start_serial, MAX(serial) AS end_serial,
        COUNT(*) AS total, SUM(CASE WHEN status == 'D' THEN 1 ELSE 0 END) AS run,
        SUM(CASE WHEN posterior != -1 THEN 1 ELSE 0 END) AS posterior
    FROM job
    GROUP BY smcSet ORDER BY smcSet;

/* smc sets where total == run AND posterior > 0 are considered complete */
/* smc sets where total == run AND posterior == 0 are considered TODO */
/* smc sets where total > run => still to be simulated */
/* TODO: if configuration data were available in db, could assess if the
*  "right" number of jobs were run, and if the "right" number of jobs have
*  posterior > -1
*/

/* add the "complete posterior" smc waves view */
CREATE VIEW IF NOT EXISTS post_done AS
    SELECT serial, smcSet, posterior
    FROM job
    JOIN (SELECT smcSet FROM smc_summary WHERE total == run AND posterior != 0)
    USING (smcSet)
    ORDER BY smcSet, posterior;

/* add the "last complete posterior" view */
CREATE VIEW IF NOT EXISTS post_last AS
    SELECT * FROM post_all JOIN (
        SELECT MAX(smcSet) AS smcSet FROM post_all
    ) USING (smcSet)
    ORDER BY posterior;

/* add the "still to process" view */
/* NB: not aware of uses were there are multiple done-but-not processed waves, but
*  can imagine a situation where that might be useful? Hence still introducing
* the all (post_todo) vs one (post_next) distinction
*/
CREATE VIEW IF NOT EXISTS post_todo AS
    SELECT serial, smcSet
    FROM job
    JOIN (SELECT smcSet FROM smc_summary WHERE total == run AND posterior == 0)
    USING (smcSet)
    ORDER BY smcSet, posterior;

/* add the "next to process" view */
CREATE VIEW IF NOT EXISTS post_next AS
    SELECT * FROM post_todo JOIN (
        SELECT MIN(smcSet) AS smcSet FROM post_todo
    ) USING (smcSet)
    ORDER BY posterior;