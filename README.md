# How to run:

- Build docker image with `docker build -t [image_name] .`
- Run docker image and portforward whichever port you would like to use for accessing the dashboard: 
    * `docker build -t [image_name] .`
    * `docker run -h localhost -p 8077:8050 -d --name [container_name] [image_name]`
        * The second number should always be 8050 because that is the port that dash will run the app on, inside the docker container. The first number is the number you will use to access the dashboard in your browser. For the first number, make sure to use a port that is not already in use by any other service (I use 8077 so it doesn't interfere with other dash apps I may be running).
    * Example:
```
docker build -t NCOG_dash_image .
docker run -h localhost -p 8077:8050 -d --name NCOG_dash_app NCOG_dash_image
```

Navigate to **localhost:8077** in a web browser (or whatever port you are using)