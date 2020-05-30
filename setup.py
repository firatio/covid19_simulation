import setuptools

setuptools.setup(
    name="covid19_sim",
    version="0.1",
    author="Firat Sabancioglu",
    description="Covid-19 Simulation",
    long_description='An agent-based simulation to discover which strategies are best to contain a covid-19 outbreak and minimize damage to freedoms',
    long_description_content_type="text/markdown",
    url="https://github.com/firatio/covid19-sim",
    packages=setuptools.find_packages(),
    keywords=['covid-19', 'epidemic simulation', 'simulation', 'computer experiment', 'agent-based model'],
    classifiers=[
        "Programming Language :: Python :: 3.6",
    ],
    python_requires='>=3.6.5',
)