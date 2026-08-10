"""

prototype to scrape and get scec benchmark data


"""
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
import pandas as pd
import os
import argparse


class SCEC:
    url = f"https://strike.scec.org/cvws/cgi-bin/seas.cgi"


    def __init__(self, headless=False, driver=None):
        self.benchmark_name = ""
        self.run_name = ""
        self.data_name = ""
        self.current_table_names = []
        self.data = pd.DataFrame()
        options = Options()
        if headless:
            options.add_argument("--headless=new")
        if driver is not None:
            self.driver = driver
            self.isopen = True
        else:
            self.driver = webdriver.Firefox(options=options)
            self.isopen = False

    def close(self):
        self.driver.quit()
        self.isopen = False
    
    def open(self):
        self.driver.get(self.url)
        self.isopen = True

    def go_back_bench(self):
        print("Going back to benchmarks page")
        self.driver.find_element(by="name",value="G0010").click()


    def go_back_user(self):
        print("Going back to runs page")

        self.driver.find_element(by="name",value="G1044").click()


    def go_back_file(self):
        print("Going back to files page")

        self.driver.find_element(by="name",value="G1046").click()


    def get_table_names(self):
        tables = self.driver.find_elements(By.TAG_NAME, value="table")
        n_tables = len(tables)
        self.current_table_names = []
        if n_tables > 1:
            for i in range(2, n_tables):
                for t in tables[i].text.splitlines()[1:]:
                    self.current_table_names.append(t.split()[0])
        else:
            for t in tables[0].text.splitlines()[1:]:
                    self.current_table_names.append(t.split()[0])


    def login(self, username=None, password=None):
        if self.isopen == False:
            self.open()

        if username is None or password is None:
            print("Username or password not provided, accessing public area")
            self.driver.find_element(by="name",value="G0012").click()
        else:
            print("Username and password provided, logging in")
            self.driver.find_element(by="name",value="G0011").click() # user login to access data button
            self.driver.find_element(by="name",value="u").send_keys(username)
            self.driver.find_element(by="name",value="p").send_keys(password)
            self.driver.find_element(by="name",value="blogin").click()

        self.get_table_names()



    def select_benchmark(self, benchmark_name):
        self.get_table_names()

        if benchmark_name in self.current_table_names:
            print(f"Selecting {benchmark_name}")

            self.benchmark_name = benchmark_name
            benchmark_name = f"G1045{benchmark_name.lower()}"
            self.driver.find_element(by="name",value=benchmark_name).click()
            self.get_table_names()
        else:
            print("benchmark not present in list")


    def select_run(self, run_name):
        self.get_table_names()

        if run_name in self.current_table_names:
            print(f"Selecting {run_name}")

            self.run_name = run_name
            run_name = f"G1047{run_name.lower()}"
            self.driver.find_element(by="name",value=run_name).click()
            self.get_table_names()
        else:
            print("run not present in list")


    def select_data(self, data_name):
        data_name = data_name.lower()
        self.get_table_names()

        if data_name in self.current_table_names:
            print(f"Selecting {data_name}")

            self.data_name = data_name
            match data_name:
                case "rupture":
                    button="G1063"
                case _:
                    button="G1059"

            for td in self.driver.find_elements(By.TAG_NAME, value="td"):

                if data_name in td.text:
                    data_name = f"{button}{data_name}"
                    self.driver.find_element(by="name",value=data_name).click()
                    return
        else:
            print("data not present")

    def parse_data(self):
        print(f"Parsing page")
        
        content = self.driver.find_element(By.TAG_NAME, value="pre").text.splitlines()
        header = None
        data_start_idx = 0 


        
        for i, c in enumerate(content):
            if "t" in c.split(): # Headers always have time as variable
                header = c.split()
                # self.data = pd.DataFrame(columns=header)

            try: # Using try to check if the components of the string are numbers -> start of the data stream
                test = c.split()
                float(test[0])
                data_start_idx = i 
                break
            except ValueError:
                ...
        data = [c.split() for c in content[data_start_idx:]]
        self.data = pd.DataFrame(data)
        self.data.columns = header


    def save_data(self, outputdir=""):

        path = os.path.join(outputdir, self.benchmark_name, self.run_name)

        if not os.path.exists(path):
            os.makedirs(path)

        filename = os.path.join(path, f"{self.data_name}.csv")
        print(f"Saving: {filename}")


        self.data.to_csv(filename, index=False, sep=";")






parser = argparse.ArgumentParser(description="SCEC SEAS benchmarks download utility")


parser.add_argument("--n", default=None, help="Name of the benchmark to download")
parser.add_argument("--r", default=None, help="Name of the run to download")
parser.add_argument("--d", nargs='+', default=[], help="Name of the data to download")
parser.add_argument("--a", action='store_true', help="Download all")
parser.add_argument("--s", default=".", help="Directory to save")
parser.add_argument("--up", nargs='+', default=[], help="scec username and password. If not provided env variables will be checked ('scec_user' and 'scec_pass'). If not present public area will be entered")



args = parser.parse_args()


scraper = SCEC(headless=False)

if len(args.up) > 1:
    scraper.login(args.up[0], args.up[1])

else:
    if "scec_user" in os.environ.keys() or "scec_pass" in os.environ.keys():
        scraper.login(os.environ["scec_user"], os.environ["scec_pass"])
    else:
        scraper.login()
if args.n: 
    scraper.select_benchmark(args.n)
else:
    scraper.select_benchmark("bp7-fd-a")

if args.r:
    scraper.select_run(args.r)
else:
    scraper.select_run("lambert")

if args.a:
    data_names = scraper.current_table_names
else:
    if args.d:
        data_names = args.d
    else:
        data_names = ["fltststrk+000dp+000"]

for name in data_names:
    scraper.select_data(name)
    scraper.parse_data()
    scraper.save_data(args.s)
    scraper.go_back_file()







