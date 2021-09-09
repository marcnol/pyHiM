
# packaging pyHiM

## Pacl using setup.py

```sh
python3 setup.py sdist bdist_wheel
```



## Pack using docker:


```sh
sudo docker build -t py_him .

sudo docker save py_him >dist/docker_pyHiM.tar
gzip dist/docker_pyHiM.tar
```

### Deploy docker container

copy docker_pyHiM.tar.gz into server, then run

```sh
docker import dist/docker_pyHiM.tar.gz new_py_him:latest
then run by :
docker run py_him
```

