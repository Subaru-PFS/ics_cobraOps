from setuptools import setup, find_packages

setup(name="ics.cobraOps",
      #version="x.y",
      author="Javier Gracia Carpio",
      #author_email="",
      #description="",
      url="https://github.com/Subaru-PFS/ics_cobraOps/",
      packages=find_packages("python"),
      package_dir={'':'python'},
      zip_safe=True,
      license="GPLv3",
      install_requires=["numpy", "matplotlib"],
      )
