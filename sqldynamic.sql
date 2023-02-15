
/* TODO: figure out how to do this portably? for now, need to do next if manually creating */
-- SELECT load_extension('sqdb/eval');

/* these blocks use `eval` to look at the contents of the (par|met)_name tables, and
* then construct a view of the (par|mar)_under tables (which are long format),
* in long format, with the column names.
* 
* The view is dropped / recreated every time this is run: the view is dynamic, but its construction is
* not: so if e.g. parameters have been added, the view columns are not updated.
*
* In generally, this step only needs to be run once, when the database is first created. If the configuration
* updates, then view may need to be recreated. However, it's also harmless (cycles aside) to run this step every time.
*/

DROP VIEW IF EXISTS par;
WITH expr AS (
    SELECT 'CREATE VIEW par AS ' AS line -- create the pivot view
    UNION ALL
    SELECT 'SELECT serial' -- pivot on serial vs...
    UNION ALL
    SELECT ', seed' -- add the seed column (from seeds table)
    UNION ALL -- all of the unique parameter indices, and label them by name
    SELECT ', MAX(value) FILTER (WHERE parIdx == ' || parIdx || ') AS ' || name FROM par_name
    UNION ALL
    SELECT ' FROM par_under JOIN seeds USING(serial) GROUP BY serial;' -- generate the pivot from the par_under(lying) table
)
SELECT eval(GROUP_CONCAT(line, '')) FROM expr; 

/* now we do the same for metrics: */

DROP VIEW IF EXISTS met;
WITH expr AS (
    SELECT 'CREATE VIEW met AS ' AS line -- create the pivot view
    UNION ALL
    SELECT 'SELECT serial ' -- pivot on serial vs...
    UNION ALL -- all of the unique metric indices, and label them by name
    SELECT ', MAX(value) FILTER (WHERE metIdx == ' || metIdx || ') AS ' || name || ' ' FROM met_name
    UNION ALL
    SELECT 'FROM met_under GROUP BY serial;' -- generate the pivot from the par_under(lying) table
)
SELECT eval(GROUP_CONCAT(line, '')) FROM expr;