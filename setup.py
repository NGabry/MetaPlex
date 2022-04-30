from setuptools import setup

setup(
    name='Metaplex',
    version='1.1.0',
    description='Read Processing and Quality Control Toolkit for Dual-Indexed Metabarcoding',
    url='https://github.com/NGabry/MetaPlex',
    author='Nick Gabry',
    author_email='n.t.gabry@gmail.com',
    license='BSD 3-clause',
    packages=['metaplex'],
    entry_points={
        'console_scripts': [
            'Metaplex-remultiplex=metaplex.remultiplexing:main',
            'Metaplex-calculate-IJR=metaplex.index_jump:main',
            'Metaplex-per-sample-filter=metaplex.per_sample_filtering:main',
            'Metaplex-length-filter=metaplex.length_filtering:main',
        ],
    },
    install_requires=['numpy',
                      'pandas',
                      'cutadapt',
                      'biom-format',
                      'biopython',
                      ],

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

    ],
)
