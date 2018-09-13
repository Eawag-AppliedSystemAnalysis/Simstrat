if [ $(docker inspect -f '{{.State.Running}}' simstrat) = "false" ]; then
        echo 'simstrat container not running... I wake it up'
        docker start simstrat
fi

docker exec simstrat sh -c "cd /home/Simstrat/Simstrat_v2/; ./run.sh"