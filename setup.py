from setuptools import setup, find_packages

setup(name="ics.cobraOps",
      version="1.0.0",
      author="Javier Gracia Carpio",
      author_email="jgracia@mpe.mpg.de",
      description="A pacakage that simulates the movements of cobras and calculates their possible collisions.",
      url="https://github.com/Subaru-PFS/ics_cobraOps/",
      packages=find_packages("python"),
      package_dir={'':'python'},
      zip_safe=True,
      license="GPLv3",
      install_requires=["numpy", "matplotlib"],
      )
