import setuptools

setuptools.setup(
    name="simcmb",
    version="0.1",
    author="Samuel D. McDermott",
    author_email="samueldmcdermott@gmail.com",
    description="code for simulating and delensing the cmb",
    packages=["simcmb"],
    url="https://github.com/deepskies/simcmb",
    install_requires=["setuptools", "camb", "namaster"],
    package_data={
        "simcmb": ["inifiles/*"],
    }
)