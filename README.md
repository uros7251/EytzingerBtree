# Optimizing B+ tree Node Layout for Improved Search Performance Using SIMD Instructions

## Project description
### Goal of the Research
The goal of this research is to enhance the performance of B+ tree structures, commonly used in
database systems, by optimizing the layout of B+ tree nodes and leveraging Single Instruction,
Multiple Data (SIMD) instructions to accelerate search operations. By adopting an implicit $n$-ary
tree layout and aligning it with SIMD capabilities, we aim to improve the efficiency of B+ tree-based data retrieval.

### Introduction
B+ trees are vital in computer science, especially within database systems, where they efficiently
index large datasets. Their balanced hierarchical structure enables swift data access, essential for
optimizing database query performance. However, conventional B+ tree node implementations
face challenges. Using sorted arrays for key storage and employing binary search during lookup
operations leads to random access patterns, hindering efficient cache utilization. Additionally, the lack of temporal locality, with most frequently accessed keys scattered across different cache lines, poses a performance problem. Finally, modern hardware supports parallel comparisons, a feature that cannot be fully utilized in a conventional sorted array setup. These issues highlight the need for an improved node implementation to enhance B+ tree efficiency.

### Proposed Modifications and Their Implications
In this research, we explore an alternative approach to optimize the layout of B+ tree nodes. Just as B+ trees were initially devised to handle latency disparities between disk and main memory, our goal is to redesign B+ tree nodes to be cache-conscious and fully leverage the SIMD capabilities of modern processors. Our approach adopts a static, pointer-free tree layout, akin to Eytzinger layout used in common binary heap implementations. The key consideration is selecting an optimal maximum number of children, denoted as $n$, in such a way that $n-1$ (the number of keys) aligns with the maximal capacity of SIMD compare instructions supported by the target machine. This alignment would enable the execution of $n−1$ comparisons in a single SIMD instruction, potentially leading to performance gains.

Additionally, this modified layout offers another potential advantage — improved temporal
locality. By arranging nodes inside a single page based on their depth within a tree and ensuring
proper alignment, they can be grouped to fit into a single cache line. This arrangement enhances
cache coherency, potentially reducing cache misses during search operations and improving search
performance.

While this layout offers theoretical prospect of performance improvements, our research will
include comprehensive empirical measurements on modern hardware to assess the actual practical
impact of these modifications accurately.

### Practical importance
Optimizing B+ tree structures for efficient data retrieval has substantial real-world implications,
particularly in the domain of database systems. If our approach proves successful in enabling
faster search operations within nodes, it could significantly enhance query execution speed, thereby improving the overall efficiency of database operations. Furthermore, faster search might offer
an opportunity to accommodate larger page sizes, potentially leading to even better performance
outcomes. In today’s data-driven world, where databases underpin a wide range of applications,
from web services to enterprise solutions, any enhancement in database efficiency holds significant practical value.

# README Template

## Getting started

To make it easy for you to get started with GitLab, here's a list of recommended next steps.

Already a pro? Just edit this README.md and make it your own. Want to make it easy? [Use the template at the bottom](#editing-this-readme)!

## Add your files

- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://gitlab.db.in.tum.de/uros7251/guided-research.git
git branch -M main
git push -uf origin main
```

## Integrate with your tools

- [ ] [Set up project integrations](https://gitlab.db.in.tum.de/uros7251/guided-research/-/settings/integrations)

## Collaborate with your team

- [ ] [Invite team members and collaborators](https://docs.gitlab.com/ee/user/project/members/)
- [ ] [Create a new merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)
- [ ] [Automatically close issues from merge requests](https://docs.gitlab.com/ee/user/project/issues/managing_issues.html#closing-issues-automatically)
- [ ] [Enable merge request approvals](https://docs.gitlab.com/ee/user/project/merge_requests/approvals/)
- [ ] [Set auto-merge](https://docs.gitlab.com/ee/user/project/merge_requests/merge_when_pipeline_succeeds.html)

## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing(SAST)](https://docs.gitlab.com/ee/user/application_security/sast/)
- [ ] [Deploy to Kubernetes, Amazon EC2, or Amazon ECS using Auto Deploy](https://docs.gitlab.com/ee/topics/autodevops/requirements.html)
- [ ] [Use pull-based deployments for improved Kubernetes management](https://docs.gitlab.com/ee/user/clusters/agent/)
- [ ] [Set up protected environments](https://docs.gitlab.com/ee/ci/environments/protected_environments.html)

***

# Editing this README

When you're ready to make this README your own, just edit this file and use the handy template below (or feel free to structure it however you want - this is just a starting point!). Thank you to [makeareadme.com](https://www.makeareadme.com/) for this template.

## Suggestions for a good README
Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information.

## Name
Choose a self-explaining name for your project.

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
