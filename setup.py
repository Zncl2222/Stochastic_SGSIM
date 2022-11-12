import setuptools

setuptools.setup(
    name='uc_sgsim',
    version='0.1.2',
    author='Zncl2222',
    author_email='3002shinning@gmail.com',
    description='Geo-Statistic alogrithm (Sequential Gaussian Simulation)',
    url='https://github.com/Zncl2222/Stochastic_UC_SGSIM',
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={'uc_sgsim': ['uc_sgsim.dll']},
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    python_requires='>=3.7',
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
    ],
)
