import numpy as np
from sklearn.cluster import DBSCAN
from NucleusModel import ConstructNucleus
from FindCase3 import AdjacencyDatabase
from collections import defaultdict
from tqdm import tqdm

fSPointsProb = 0.17
fEMinDamage = 5.0
fEMaxDamage = 37.5
fMinDist = 3.2
minTADbp = 100000

def rotate_points(points, theta_rad, phi_rad):
    # 定义绕z轴旋转的旋转矩阵
    Rz = np.array([
        [np.cos(theta_rad), -np.sin(theta_rad), 0],
        [np.sin(theta_rad), np.cos(theta_rad), 0],
        [0, 0, 1]
    ])

    # 定义绕新的y轴旋转的旋转矩阵
    Ry = np.array([
        [np.cos(phi_rad), 0, np.sin(phi_rad)],
        [0, 1, 0],
        [-np.sin(phi_rad), 0, np.cos(phi_rad)]
    ])

    # 先绕z轴旋转，再绕新的y轴旋转
    rotated_points = points.dot(Rz).dot(Ry)

    return rotated_points

class HitPoints:
    def __init__(self):
        self.DNAModel = ConstructNucleus()

    def AddOnePoints(self, x,y,z,E,id):
        self.id_read.append(id)
        self.position_read.append([x,y,z])
        self.energy_read.append(E)

    def BeforeReadAll(self):
        self.id_read = []
        self.position_read = []
        self.energy_read = []

    def ReadHitPoints(self, filename, bias=0):
        with open(filename, 'r') as f:
            with open(filename, 'r') as file:
                for line in file:
                    # 忽略空行
                    if not line.strip():
                        continue
                    # 假设每行的数据格式为 id E x y z
                    parts = line.strip().split()
                    id = int(parts[0]) + bias
                    E = float(parts[1])  # 假设E是浮点数
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])  # 假设x, y, z都是浮点数
                    self.AddOnePoints(x,y,z,E,id)

    def AfterReadAll(self):
        self.position_read = np.array(self.position_read)
        self.energy_read = np.array(self.energy_read)
        self.id_read = np.array(self.id_read)

    def AssignStrand(self):
        rand = np.random.uniform(0, 1.0, size=len(self.energy))
        self.strand = np.where(rand > 0.5, 1, 2)

    def IsEdepSufficient(self):
        proba = (self.energy - fEMinDamage) / (fEMaxDamage - fEMinDamage)
        proba[self.energy<fEMinDamage] = 0.0
        proba[self.energy>fEMaxDamage] = 1.0
        rand = np.random.uniform(0, 1.0, size=len(proba))
        result = proba > rand
        return result

    def DeleteHitsNotInduceDamage(self):
        EdepSufficient = self.IsEdepSufficient()
        self.position = self.position[EdepSufficient]
        self.id = self.id[EdepSufficient]
        self.energy = self.energy[EdepSufficient]
        rand = np.random.uniform(0, 1.0, size=len(self.energy))
        HitDNA = rand < fSPointsProb
        mask = HitDNA
        self.position = self.position[mask]
        self.id = self.id[mask]
        self.energy = self.energy[mask]
        InTAD3 = self.DNAModel.IsPointsInTAD(self.position)
        mask = InTAD3
        self.position = self.position[mask]
        self.id = self.id[mask]
        self.energy = self.energy[mask]
        self.AssignStrand()


    def FindClusters(self):
        dbscan = DBSCAN(eps=fMinDist, min_samples=1, metric='euclidean')
        dbscan.fit(self.position)
        labels = dbscan.labels_
        return labels

    def FindDSBs(self, labels):
        self.dsb_clusters = []
        ssb_clusters = []
        dsb_count = 0
        ssb_count = 0
        for cluster_id in set(labels):
            if cluster_id == -1:
                # -1是噪声点的标签，忽略它
                continue

            # 找出属于当前cluster的所有点的索引
            cluster_indices = np.where(labels == cluster_id)[0]
            # 检查当前cluster中的strand值
            cluster_strands = self.strand[cluster_indices]
            # 如果cluster中同时包含1和2，则为DSB
            if 1 in cluster_strands and 2 in cluster_strands:
                self.dsb_clusters.append(cluster_indices)
                dsb_count += 1
            else:
                ssb_clusters.append(cluster_indices)
                ssb_count += 1
        return dsb_count, ssb_count

    def AnalyzeDSB(self):
        particles = []
        TADs = []
        for indices in self.dsb_clusters:
            ids = self.id[indices]
            pos = self.position[indices]
            TAD = self.DNAModel.InWhichTAD(pos[0])
            if self.DNAModel.TADs[TAD].bp > minTADbp:
                TADs.append(TAD)
                particles.append(ids[0])

        # Step 1: Group TADs into clusters
        clusters = defaultdict(list)
        for index, tad in enumerate(TADs):
            clusters[tad].append(index)

        cluster_used = clusters
        # 首先寻找case3，有频繁作用的TAD就算做一个case3
        db = AdjacencyDatabase()
        for TAD in cluster_used.keys():
            for neighbor in self.DNAModel.TADs[TAD].neighbors:
                db.add(TAD, neighbor)  # 确保邻居也被添加到并查集中

        case3_cluster = db.find_case3_cluster(cluster_used.keys())

        count_TAD2 = 0
        count_TAD2_dif = 0
        for cluster in case3_cluster:
            if len(cluster) >= 2:
                count_TAD2 += 1
                unique_particles = set()
                for hit_p in cluster:
                    indices = cluster_used[hit_p]
                    for index in indices:
                        unique_particles.add(particles[index])
                if len(unique_particles) > 1:
                    count_TAD2_dif += 1

        for cluster in case3_cluster:
            if len(cluster) >= 2:
                for hit_p in cluster:
                    cluster_used.pop(hit_p, None)

        # Step 2: Count clusters with more than one element
        clusters_with_more_than_one_element = sum(1 for indices in cluster_used.values() if len(indices) > 1)

        # Step 3: Count clusters with different particles
        clusters_with_different_particles = 0
        for indices in cluster_used.values():
            if len(indices) > 1:
                unique_particles = set(particles[index] for index in indices)
                if len(unique_particles) > 1:
                    clusters_with_different_particles += 1

        cluster_one = {}
        for TAD, indices in cluster_used.items():
            if len(indices) == 1:
                cluster_one[TAD] = indices



        return len(cluster_used), clusters_with_more_than_one_element, clusters_with_different_particles, count_TAD2, count_TAD2_dif

    def RunOneTime(self, filter=1000000):
        self.id = self.id_read[self.id_read<filter]
        self.position = self.position_read[self.id_read<filter]
        self.energy = self.energy_read[self.id_read<filter]

        theta = np.arccos(np.random.uniform(-1, 1))
        phi = np.random.uniform(0.0, 2.0 * np.pi)
        self.position = rotate_points(self.position, theta, phi)

        self.DeleteHitsNotInduceDamage()
        labels = self.FindClusters()
        dsb, ssb = self.FindDSBs(labels)
        cluster, cluster_plus, cluster_plus_two, hit_2TAD, hit_2TAD_dif = self.AnalyzeDSB()
        all_result = [dsb, ssb, cluster, cluster_plus, cluster_plus_two, hit_2TAD, hit_2TAD_dif]
        return all_result



if __name__ == "__main__":
    hits = HitPoints()

    for iter in range(0,20):
        print(f"iter {iter}")
        hits.BeforeReadAll()
        for run in range(iter*10, (iter+1)*10):
            hits.ReadHitPoints(f"electron/{run}/Trackstructure_ID.txt", bias=5000 * (run-iter*10))
            hits.ReadHitPoints(f"electron/{run}/Subexcitation_ID.txt", bias=5000 * (run-iter*10))
        hits.AfterReadAll()

        run_times = 10
        for particle in range(0, 10):
            damage_result = np.zeros([run_times, 7], dtype=int)
            for i in tqdm(range(run_times)):
                aResult = hits.RunOneTime(filter=(particle + 1) * 5000)
                damage_result[i, :] = aResult

            output_name = f'electron_run017-2/dsb_result{iter}_{particle}.txt'
            with open(output_name, 'ab') as file:
                np.savetxt(file, damage_result, fmt='%d', delimiter='\t')
