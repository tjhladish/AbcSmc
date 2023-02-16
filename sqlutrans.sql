
/* normal-formed untransformed parameters table. links run serial + parameter id => untrans value of that parameter */
CREATE TABLE IF NOT EXISTS upar_vals (
    serial INTEGER NOT NULL,
    parIdx INTEGER NOT NULL,
    value REAL,
    PRIMARY KEY (serial, parIdx)
);

/* mirror of par view, but for untransformed parameters */
DROP VIEW IF EXISTS upar;
WITH expr AS (
    SELECT 'CREATE VIEW upar AS ' AS line -- create the pivot view
    UNION ALL
    SELECT 'SELECT serial' -- pivot on serial vs...
    UNION ALL
    SELECT ', COALESCE(seed, serial) AS seed' -- add a seed column (from seeds table if filled, but defaults to serial)
    UNION ALL -- all of the unique parameter indices, and label them by name
    SELECT ', MAX(value) FILTER (WHERE parIdx == ' || parIdx || ') AS ' || name FROM par_name
    UNION ALL
    SELECT ' FROM upar_vals LEFT JOIN seeds USING(serial) GROUP BY serial;' -- generate the pivot from the par_vals(lying) table
)
SELECT eval(GROUP_CONCAT(line, '')) FROM expr; 
