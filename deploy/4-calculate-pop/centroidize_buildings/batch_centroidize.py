#centroidize_deploy.py
import argparse
from centroidize import centroidize

def setup():
	parser = argparse.ArgumentParser(description='Construct building centroids and areas from building files')
	parser.add_argument('--top_building_dir', dest='top_building_dir', type=str, required=True)
	return parser.parse_args()

if __name__ == '__main__':
	centroidize(**vars(setup()))