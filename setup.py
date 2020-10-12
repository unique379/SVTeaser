from setuptools import setup

setup(
    name='SVTeaser',
    version='0.0.1',
    author="ACEnglish",
    author_email="acenglish@gmail.com",
    url="https://github.com/collaborativebioinformatics/SVTeaser",
    packages=['svteaser'],
    license='MIT',
    scripts=["bin/svteaser"],
    description="SV simulation for rapid benchmarking",
    long_description=open('README.md', encoding='UTF-8').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        "truvari>=2.0.1"
    ],
)
