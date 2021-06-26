ANALYZE;

DELETE
FROM `sequence`
WHERE `ncbi_id` NOT IN
    (SELECT `ncbi_id`
     FROM `taxonomy`
     WHERE `name_class` = 'scientific name'
       AND `node_rank` = 'species'
       AND `left_value` > (SELECT `left_value` FROM `taxonomy` WHERE `name` = 'Actinopterygii')
       AND `right_value` < (SELECT `right_value` FROM `taxonomy` WHERE `name` = 'Actinopterygii')
    );
-- QUERY PLAN
-- |--SCAN TABLE SEQUENCE
-- `--LIST SUBQUERY 3
--    |--SEARCH TABLE taxonomy USING INDEX taxonomy_right_value (right_value<?)
--    |--SCALAR SUBQUERY 2
--    |  `--SEARCH TABLE taxonomy USING INDEX taxonomy_name (name=?)
--    `--SCALAR SUBQUERY 1
--       `--SEARCH TABLE taxonomy USING INDEX taxonomy_name (name=?)

DELETE
FROM `taxonomy`
WHERE `left_value` < (SELECT `left_value` FROM `taxonomy` WHERE `name` = 'Actinopterygii')
   OR `right_value` > (SELECT `right_value` FROM `taxonomy` WHERE `name` = 'Actinopterygii');

-- QUERY PLAN
-- `--MULTI-INDEX OR
--    |--INDEX 1
--    |  |--SCALAR SUBQUERY 1
--    |  |  `--SEARCH TABLE taxonomy USING INDEX taxonomy_name (name=?)
--    |  `--SEARCH TABLE taxonomy USING INDEX taxonomy_left_value (left_value<?)
--    `--INDEX 2
--       |--SCALAR SUBQUERY 2
--       |  `--SEARCH TABLE taxonomy USING INDEX taxonomy_name (name=?)
--       `--SEARCH TABLE taxonomy USING INDEX taxonomy_right_value (right_value>?)

VACUUM;
ANALYZE;
VACUUM;
