class Variants:
    def __init__(self, project, rna=False, filtered=False):
        pass

    def display_impacts(self):
        pass

    def filter(self, filters):
        pass

    def fields(self):
        pass

    def __str__(self):
        pass




class Variant:
    """
    same class will be used for both rna and dna variants
    """
    def __init__(self, project):
        pass

    def calc_freq(self, samples=None, cohort=None):
        pass

    def find_samples(self, cohort=None, genotype="both"):
        pass

    def impact(self, highest=True):
        pass

    def __str__(self):
        pass