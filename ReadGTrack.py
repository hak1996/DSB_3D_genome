
import sys
import pandas as pd

class GtrackFile(object):

    def __init__(self):
        self.header = None
        self.chromosomes = {}
        self.bead_id_chr = {}
        self.periphery_ids = []
        self.center_ids = []
        self.interaction_ids = []
        self.interaction_distance_ids = []
        self.interaction_lower_distance_ids = []
        self.interaction_upper_distance_ids = []
        self.dist_info = {}
        self.weight_info = {}
        self.periphery_weight = {}
        self.center_weight = {}
        self.bead_radius = 0.5  # default value
        self.has_periphery = False
        self.has_edges = False
        self.has_colour = False
        self.has_bead_radius = False
        self.max_radius = 0

    def read_gtrack(self, fp: str):
        f = open(fp, "r")
        data_set = {}
        for i, line in enumerate(f.readlines()):
            if line.startswith("###"):
                self.header = line.replace('###', '').replace('\n', '').split('\t')
                self.has_periphery = "periphery" in self.header
                self.has_edges = "edges" in self.header
                self.has_colour = "color" in self.header or "colour" in self.header
                self.has_bead_radius = "radius" in self.header
                data_set = {head: [] for head in self.header}
                continue
            elif line.startswith("#"):  # comments
                continue

            for header, data in zip(self.header, line.replace('###', '').replace('\n', '').split('\t')):
                data_set[header].append(data)

            # filling in this information here to save iterating through it later
            # periphery
            if self.has_periphery and data_set['periphery'][-1] != '.':
                if data_set['periphery'][-1] == '1':
                    self.periphery_ids.append(data_set['id'][-1])
                elif data_set['periphery'][-1] == '0':
                    self.center_ids.append(data_set['id'][-1])
                elif ',' in data_set['periphery'][-1]:  # has periphery weighting
                    if data_set['periphery'].split(',')[0] == '1':
                        self.periphery_ids.append(data_set['id'][-1])
                        self.periphery_weight[self.periphery_ids[-1]] = float(data_set['periphery'].split(',')[1])
                    elif data_set['periphery'].split(',')[0] == '0':
                        self.center_ids.append(data_set['id'][-1])
                        self.center_weight[self.center_ids[-1]] = float(data_set['periphery'].split(',')[1])
                    else:
                        sys.exit(f'Gtrack Line Error {i}: Periphery can be either 0 or 1, when a weight is given.')

            # edges
            if self.has_edges and data_set['edges'][-1] != '.':
                edge_list = data_set['edges'][-1].split(';')
                for edge in edge_list:
                    if edge == data_set['id'][-1] or edge.split('=')[0] == data_set['id'][-1]:
                        sys.exit(f'Gtrack Line Error {i}: Same edge value as own BeadID.')
                    elif '=' not in edge:
                        id_pair = (data_set['id'][-1], edge)
                        if id_pair not in self.interaction_ids:
                            self.interaction_ids.append(id_pair)
                    else:
                        edge_id = edge.split('=')[0]
                        edge_weight = edge.split('=')[1]
                        id_pair = (data_set['id'][-1], edge_id)
                        if '.' not in edge_weight:  # No weighting actually applied
                            self.interaction_ids.append(id_pair)
                        else:
                            edge_details = edge_weight.split(',')
                            if edge_details[2] not in ["0", "1", "."]:
                                sys.exit(f'Gtrack Line Error {i}: Boundary information value can only be 0, 1 or .')
                            elif float(edge_details[0]) < 0 or float(edge_details[1]) < 0:
                                sys.exit(f'Gtrack Line Error {i}: Weight and Distance should be greater than 0')
                            elif edge_details[2] == '0':
                                self.interaction_upper_distance_ids.append(id_pair)
                            elif edge_details[2] == '1':
                                self.interaction_lower_distance_ids.append(id_pair)
                            elif edge_details[2] == '.':
                                self.interaction_distance_ids.append(id_pair)
                            self.weight_info[id_pair] = float(edge_details[0])
                            self.dist_info[id_pair] = float(edge_details[1])

        return self.interaction_ids


if __name__ == "__main__":
    #gtrack_file = 'IMR90_noLADs.gtrack'
    #gtrack = GtrackFile()
    #interaction_ids = gtrack.read_gtrack(gtrack_file)

    file_dir = "IMR90_noLADs_1.parquet"
    df = pd.read_parquet(file_dir, engine='pyarrow')

    beadID = df.beadID.tolist()
    bp = df.bp.tolist()
    x = df.x.tolist()
    y = df.y.tolist()
    z = df.z.tolist()
    r = df.radius.tolist()
    radius = [float(item) for item in r]
    print("success!")
