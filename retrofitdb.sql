
/* add the "smc set" view */
CREATE VIEW IF NOT EXISTS smc_summary AS
    SELECT smcSet, MIN(serial) AS start_serial, MAX(serial) AS end_serial,
        COUNT(*) AS total, SUM(CASE WHEN status == 'D' THEN 1 ELSE 0 END) AS run,
        SUM(CASE WHEN posterior != -1 THEN 1 ELSE 0 END) AS posterior
    FROM job
    GROUP BY smcSet ORDER BY smcSet;

/* smc sets where total == run AND posterior > 0 are considered complete */
/* TODO: if configuration data were available in db, could assess if the
*  "right" number of jobs were run, and if the "right" number of jobs have
*  posterior > -1
*/

/* add the "complete posterior" smc waves view */
CREATE VIEW IF NOT EXISTS post_done AS
    SELECT serial, smcSet, posterior
    FROM job WHERE posterior != -1 AND status == 'D'
    ORDER BY smcSet, posterior;

/* add the "last complete posterior" view */
CREATE VIEW IF NOT EXISTS post_last AS
    SELECT * FROM post_all JOIN (
        SELECT MAX(smcSet) AS smcSet FROM post_all
    ) USING (smcSet)
    ORDER BY posterior;

CREATE VIEW IF NOT EXISTS post_todo AS
    SELECT * FROM post_all JOIN (
        SELECT MAX(smcSet) AS smcSet FROM post_all
    ) USING (smcSet)
    ORDER BY posterior;
