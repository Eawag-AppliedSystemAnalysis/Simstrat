## How to setup Simstrat building environment with Docker
### Download the repositories
Clone the repository:

~~~bash
git clone https://github.com/Eawag-AppliedSystemAnalysis/Simstrat.git
~~~

NB: currently this feature is in the development branch `simstrat2beta`, so (before the release) jump into the branch with:

~~~bash
git checkout simstrat2beta
~~~

### Setup the Docker container
1. Install the [_docker_](https://www.docker.com/
) engine, as explained for [MacOS](https://docs.docker.com/docker-for-mac/), [Ubuntu](https://docs.docker.com/install/linux/docker-ce/ubuntu/), or [Windows](https://docs.docker.com/docker-for-windows/) users.
2. I prepared a ready-to-use *Dockerfile* with all the needed packages, so open the terminal (or run the docker toolbox) and navigate to the folder `Simstrat/Documentation/developer/docker_build_env` where the Dockerfile (and this README) is located.
3. WATCH! For Linux and (maybe) MacOS users: by default _docker engine_ requires root privileges, i.e. you have to type `sudo` in front of each of the following commands. If you wish to avoid this, you can follow this [guidelines](https://docs.docker.com/install/linux/linux-postinstall/#manage-docker-as-a-non-root-user).
4. Build the docker image using the given Dockerfile

    ~~~bash
    docker docker build -t simstrat:alpine .
    ~~~

5. Create the actual docker container (that makes use of the image you just created), remebering to provide the `<pathToLocalGitRepoDirectory>`:

    ~~~bash
    docker create --name simstrat -it -v <pathToLocalGitRepoDirectory>:/home/Simstrat simstrat:alpine
    ~~~

6. A container named `simstrat` has just been created, now let's start it:

    ~~~bash
    docker start simstrat
    ~~~

7. Your container is now running, to jump into it, open a `bash` terminal with:

    ~~~bash
    docker exec -it simstrat bash
    ~~~

8. If everything went smooth, you're now inside your virtual build environment at the path `/home/Simstrat`, which is mapping your local directory at `<pathToLocalGitRepoDirectory>`, as example.

9. Next time you want to access the container, simply execute points 6 and 7 (even just 7 if your container is already running). Now you are ready to actually build Simstrat.

### Build and run Simstrat inside the docker container
Jump into the correct branch folder and simply execute the shell script:
-  `./build_linux.sh` to build the program,
-  `./run.sh` to run the Zurichsee testcase.

### Build and run Simstrat from your hosting system:
Jump into the correct branch folder and simply execute the shell/batch script:
-  `./build_docker.sh` (or `./build_docker_Win.bat`) to build the program from a Linux (Win) host,
-  `./run_docker.sh` (or `./run_docker_Win.bat`) to run the Zurichsee testcase.

Enjoy your GLM simulations!
