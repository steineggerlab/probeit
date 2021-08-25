from setuptools import setup, find_packages

setup(
    name = 'probeit',
    version = 'v1.4',
    packages = find_packages(),
    entry_points = {
        'console_scripts': [
            'probeit = probeit.__main__:main'
        ]
    })