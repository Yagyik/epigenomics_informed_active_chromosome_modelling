from setuptools import setup, find_packages

setup(
    name='ChromUtils',
    version='0.1',
    packages=find_packages(include=['ChromUtils', 'ChromUtils.*']),
    install_requires=[
        # Add any dependencies your package needs here
    ],
    entry_points={
        'console_scripts': [
            # Add any command line scripts here
        ],
    },
    author='Your Name',
    author_email='your.email@example.com',
    description='A package for chromosome modelling utilities',
    url='https://github.com/yourusername/ChromUtils',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)