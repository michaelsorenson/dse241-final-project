# Michael's Marine Microbial Dashboard (Allen Lab)
## Prerequisites
- python libraries in requirements.txt, to install run `pip install -r requirements.txt`
- Data files for 16S, 18Sv4, 18Sv9, and metadata
    - 16S: `NCOG_21_16S_redo2_asv_count_tax.tsv`
    - 18Sv4: `NCOG_18sV4_asv_count_tax.tsv`
    - 18Sv9: `NCOG_18sV9_asv_count_tax_S.tsv`
    - Metadata: `NCOG_sample_log_DNA_stvx_meta_2014-2020_mod`

## Usage
1. Install the required libraries using `pip install -r requirements.txt`
2. Run the dashboard with `python app.py [port]`
    - Port is an optional argument, it should be a port number available on your computer, if there is no port it will default to 8050
3. Open the dashboard in your preferred web browser, `localhost:8050` or `localhost:[port]` if you have chosen a different port


# Using Docker, not complete yet:
## How to run with docker:

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