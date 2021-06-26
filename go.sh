brew install jonchang/biology/phlawd_db_maker
phlawd_db_maker vrt vrt.db
cat reduce_db.sql | sqlite3 vrt.db
ls -lha
