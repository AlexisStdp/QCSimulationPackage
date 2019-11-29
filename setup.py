import setuptools

with open("README.md", "r") as ld:
    long_description = ld.read()

setuptools.setup(
    name="QCSimulationPackage",
    version="0.0.1",
    author="AlexisStdp",
    author_email="author@example.com",
    description="Toy project simulation of a Quantum Circuit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AlexisStdp/QCSimulationPackage",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
