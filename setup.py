import setuptools
setuptools.setup(
    name="pygamma_n",
    version="0.0.1",
    author="Garrett Finucane",
    author_email="garrettdreyfus@gmail.com",
    description="A python adaptation of the MatLab Gamma N package for neutral density",
    url="https://github.com/garrettdreyfus/py-gamma-n",
    license = "MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'numpy',
        'xarray',
    ],
    packages=['pygamma_n'],
)
