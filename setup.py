from setuptools import setup, find_packages

setup(
    name='clinpy',
    version='0.1.0',
    description="Python module for creating and querying a database that has clinical -omics data",
    author='Alper Celik',
    author_email='alper.celik@sickkids.ca',
    packages=find_packages(),
    install_requires=["pandas",
                      "pyranges",
                      "sqlalchemy",
                      "numpy",
                      "pysam",
                      "pytxdb @ git+https://github.com/celalp/pytxdb@master",
                      "pyyaml"],
    zip_safe=False,
    scripts=["clinpy/scripts/create_project.py"],
    include_package_data=True
)
