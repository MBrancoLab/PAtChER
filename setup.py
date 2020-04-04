import setuptools

__pkg_name__ = "PAtChER"

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name=__pkg_name__,
    version="0.0.1",
    author="Example Author",
    author_email="author@example.com",
    description="PAtChER is a tool to help re-assign non uniquely mapping reads within a HiChIP experiment.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MBrancoLab/PAtChER",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "patcher = PAtChER.patcher:main".format(__pkg_name__),
        ],
    },
    python_requires='>=3.6',
)