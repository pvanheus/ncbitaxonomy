CREATE TABLE taxonomy (
    id INTEGER PRIMARY KEY,
    ancestry TEXT,
    name TEXT NOT NULL UNIQUE,
    rank TEXT
);

CREATE UNIQUE INDEX taxonomy_name_idx ON taxonomy(name);
