from setuptools import setup, find_packages

setup(
    name="simcmb",
    version="0.1",
    author="Samuel D. McDermott",
    author_email="samueldmcdermott@gmail.com",
    description="code for simulating and delensing the cmb",
    packages=find_packages(),
    url="https://github.com/deepskies/simcmb",
    install_requires=["setuptools", "camb", "pymaster", "astropy"],
    package_data={
        "simcmb": ["inifiles/*"],
    }
)