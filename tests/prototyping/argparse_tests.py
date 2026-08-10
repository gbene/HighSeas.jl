import argparse


parser = argparse.ArgumentParser(description="Argparse tests")

parser.add_argument("--a", default=1, help="This is the 'a' option")
parser.add_argument("--b", '--names-list', nargs='+', default=[])
parser.add_argument("--c", action='store_true', help="This is the 'c' option")


args = parser.parse_args()

print(args)