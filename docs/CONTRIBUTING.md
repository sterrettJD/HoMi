# Contributing

Welcome to the HoMi! We appreciate your interest in contributing. Before you get started, please take a moment to read and follow these contribution guidelines to ensure a smooth and collaborative development process.

## Table of Contents

1. [Getting Started](#getting-started)
    - [Fork the Repository](#fork-the-repository)
    - [Clone Your Fork](#clone-your-fork)
    - [Install HoMi](#install-dependencies)
2. [Contribution Workflow](#contribution-workflow)
    - [Branching](#branching)
    - [Coding Standards](#coding-standards)
    - [Commit Messages](#commit-messages)
3. [Testing](#testing)
4. [Submitting Pull Requests](#submitting-pull-requests)
5. [Code Reviews](#code-reviews)
6. [Documentation](#documentation)
7. [Code of Conduct](#code-of-conduct)
8. [Acknowledgments](#acknowledgments)

# Getting Started

## Fork the Repository

Click the "Fork" button on the GitHub repository to create your copy.

## Clone Your Fork

Clone the repository to your local machine:

```bash
git clone git@github.com:sterrettJD/HoMi.git
cd HoMi
```

## Install HoMi
```bash
pip install -e .
```

# Contribution Workflow
## Branching
Create a new branch for your contribution:

```bash
git checkout -b new-feature
```

## Coding Standards
Please follow PEP8 standards for Python, and use informative variable, function, and rule names. Add docstrings for new rules. Please feel free to reach out to John Sterrett for style questions. The style guide for this project is currently in development.

## Commit Messages
Write clear and concise commit messages. Please reference any relevant issues or pull requests in your commits.

# Testing
Ensure that your changes do not break existing functionality. Write tests for new features and make sure all tests pass before submitting a pull request.

```bash
python -m pytest
```

# Submitting Pull Requests
When you're ready, push your changes to your forked repository and submit a pull request to the main branch of the original repository. Provide a detailed description of your changes and reference any related issues.

# Code Reviews
All contributions go through a code review process. Please be responsive to feedback and be prepared to make changes to your code if necessary.

# Documentation
Update the documentation if your changes affect the usage or configuration of the HoMi pipeline. Documentation is as important as code!

# Code of Conduct

## Our Standards

In the interest of fostering an open and welcoming community, we as contributors and maintainers pledge to make participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.


Examples of behavior that contributes to creating a positive environment include:

- Being respectful of differing viewpoints and experiences
- Gracefully accepting constructive criticism
- Focusing on what is best for the community
- Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

- The use of sexualized language or imagery and unwelcome sexual attention or advances
- Trolling, insulting/derogatory comments, and personal or political attacks
- Public or private harassment
- Publishing others' private information, such as a physical or electronic address, without explicit permission
- Other conduct that could reasonably be considered inappropriate in a professional setting

## Responsibilities
Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

## Scope
This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project email address, posting via an official social media account, or acting as an appointed representative at an online or offline event.

## Enforcement
Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting John Sterrett (sterrettJD @ protonmail dot com). All complaints will be reviewed and investigated and will result in a response that is deemed necessary and appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution
This Code of Conduct is adapted from the Contributor Covenant, version 2.0, available at https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.


### Acknowledgments
Thank you for contributing to HoMi! Your time and effort make this project better for everyone.