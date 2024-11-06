import os

import setuptools

# since we have collaberation within this project
# api should not change arbitrarily

requirements_path = os.path.join(os.path.dirname(__file__), "requirements.txt")
if os.path.exists(requirements_path):
    with open(requirements_path, "r") as f:
        install_requires = list(map(lambda x: x.strip(), f.readlines()))
        extras_require = {}
    #     extras_require["complete"] = sorted(
    #         {v
    #          for req in extras_require.values() for v in req})
else:
    install_requires = []
    extras_require = {}


def main():
    setuptools.setup(
        name="NeoTopology",
        description="NeoBinder ltd(域新说) topology toolkit",
        packages=setuptools.find_packages(),
        python_requires=">=3.7",
        extras_require=extras_require,
        install_requires=install_requires,
        scripts=["bin/fixpdb.py"],
    )


if __name__ == "__main__":
    main()
