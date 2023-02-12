#!/bin/sh

sqlite3 $1 "CREATE VIEW IF NOT EXISTS smc_view AS SELECT serial, smcSet, posterior FROM job WHERE posterior != -1 AND status == 'D' ORDER BY smcSet, posterior;"
sqlite3 $1 "CREATE VIEW IF NOT EXISTS smc_last AS SELECT * FROM smc_view JOIN (SELECT MAX(smcSet) AS smcSet FROM smc_view) USING (smcSet) ORDER BY posterior;"