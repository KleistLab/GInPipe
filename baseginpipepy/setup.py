import setuptools

VERSION = '0.0.1'
DESCRIPTION = 'GInPipe'
LONG_DESCRIPTION = 'GInPipe - set of tools used by the pipeline'

setuptools.setup(
        name="ginpipepy",
        version=VERSION,
        author="Yannick Duport, Maureen Smith, Maria Trofimova",
        author_email="<>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=setuptools.find_packages(),
        install_requires=['pysam','numpy<1.20','Bio','scipy','pandas'],
        keywords=['python'],
        classifiers= [
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)
