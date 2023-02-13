# Major refactoring branch

This branch will have some breaking changes, the order of dev milestones are as follows:

V1 release will have the following:
1. Finish variant methods:
   + This includes variant searching with arbitrary filters
   + nested filters should be allowed
   + Testing of the methods
2. Re factoring of assays into their own folder
   + support for multiple databases during creation
   + re factor the samples and subject and cohorts tables
   + There are a couple of mandatory rules:
     + A subject can belong to multiple cohorts
     + A sample cannot come from multiple subjects
   + Each assay will have 3 files (at least) the files are 
   + tables.py (this has the table definitions, this can be 
   a sqlachemy class or a dict to be converted into a table)
   + methods.py this has the query methods, usually one class or two
   + import.py these are the import functions, they will be 
   imported in the create_project.py and will be called
   based on project description
   + Project class is the main connection storage
   + The database needs to store the assays that's in it
   + The stored assay types then need to be dynamically imported
   in the project init section
   + abililty to add samples after create, there is some that's been
   done on that front but never tested

V1.5 will need to have

1. Data versioning
2. QC metrics as an assay feature like QC.py or somesuch
    + this can be a part of the import section as well
3. unit tests

V2 will need to have analysis sections and pipeline information as 
well as analysis versioning. 


V3 will need support for very large arrays (SC ATAC or RNA-Seq style)