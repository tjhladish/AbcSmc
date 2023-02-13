#!/bin/sh

sqlite3 $1 "CREATE VIEW IF NOT EXISTS post_all AS SELECT serial, smcSet, posterior FROM job WHERE posterior != -1 AND status == 'D' ORDER BY smcSet, posterior;"
sqlite3 $1 "CREATE VIEW IF NOT EXISTS post_last AS SELECT * FROM post_all JOIN (SELECT MAX(smcSet) AS smcSet FROM post_all) USING (smcSet) ORDER BY posterior;"
sqlite3 $1 "CREATE VIEW IF NOT EXISTS unfiltered ..."