import warnings

class Assay:
    def __init__(self, name, project, create=False):
        """
        table_dict is tablename:table so for example rna_variants:variants key is the name of the table value is the
        kind of table that is defined in import_data methods
        """
        self.project=project
        self.name=name
        self.create=create

    def check_tables(self, table_dict, quiet=True):
        self.project.metadata.reflect()
        current_tables=table_dict.keys()
        not_in=[]
        for table in self.tables.keys():
            if table in current_tables:
                continue
            else:
                not_in=[table]
                if not quiet:
                    warnings.warn("{} table which is a {} table is missing from the database")

        if len(not_in)==len(self.tables.keys):
            self.create=True
        elif not_in==0:
            self.create=False
        else:
            warnings.warn("{} are missing tables from {} may be there "
                          "was an error during import".format(",".join(not_in), self.name))

    def import_assay_data(self, import_fun, file, read_fun, assay_params, read_fun_params):
        import_fun(file, self.project, read_fun, assay_params, read_fun_params, create=self.create)
