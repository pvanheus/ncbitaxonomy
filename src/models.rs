use super::schema::taxonomy;

#[derive(Queryable)]
pub struct Taxon {
    pub id: i32,
    pub ancestry: Option<String>,
    pub name: String,
    pub rank: Option<String>
}

#[derive(Insertable)]
#[table_name="taxonomy"]
pub struct NewTaxon<'a> {
    pub id: &'a i32,
    pub ancestry: Option<&'a str>,
    pub name: &'a str,
    pub rank: Option<&'a str>
}