#!/bin/bash

{
    usr_id=$(sudo -u ${SUDO_USER-$USER} id -u)
} || {
    usr_id=$(id -u)
}

echo "Creating a user 'developer' with the permissions of user $usr_id."
echo

docker build --build-arg developer_user_id=$usr_id -t simstrat:ubuntu "$(pwd)/$(dirname $0)"
docker create --name simstrat -it -v "$(pwd)/$(dirname $0)/../../":/home/Simstrat simstrat:ubuntu