# -*- coding: utf-8 -*-
import torch
print(torch.__version__)
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Variable
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph, PSObject
from ElsevierAPI.ResnetAPI.ResnetAPIcache import APIcache
from torch.autograd import Variable

# Install required packages.
import os
os.environ['TORCH'] = torch.__version__


class SkipGram(nn.Module):
    def __init__(self, num_nodes, emb_size, num_features):
        super(SkipGram, self).__init__()
        self.node_embeddings = nn.Embedding(num_nodes, emb_size)
        self.feature_embeddings = nn.Embedding(num_features, emb_size)

    def forward(self, target_node, context_node, target_feature, context_feature):
        emb_target_node = self.node_embeddings.forward(target_node)
        emb_context_node = self.node_embeddings.forward(context_node)

        target_feature_tensor = torch.tensor(target_feature)
        #target_feature_tensor_clone = target_feature_tensor.clone().detach().requires_grad_(True)
        target_feature_index = torch.nonzero(target_feature_tensor)[0, 0]

        context_feature_tensor = torch.tensor(context_feature).clone().detach()
        context_feature_index = torch.nonzero(context_feature_tensor)[0, 0]

        emb_target_feature = self.feature_embeddings(torch.LongTensor([target_feature_index]))
        emb_context_feature = self.feature_embeddings(torch.LongTensor([context_feature_index]))

        scores_node = torch.matmul(emb_target_node, emb_context_node.t())
        scores_feature = torch.matmul(emb_target_feature, emb_context_feature.t())

        return scores_node + scores_feature


def load_sample_resnet(cache_name:str):
    sample = APIcache(cache_name=cache_name,connect2server=False)
    return sample.network


def random_walk(edge_list, start_node, walk_length):
    current_node = start_node
    walk = [current_node]
    
    for _ in range(walk_length):
        neighbors = [(v, w, f) for (u, v, w, f) in edge_list if u == current_node]
        if not neighbors:
            break
        probabilities = torch.tensor([w for (_, w, _) in neighbors], dtype=torch.float32)
        probabilities /= probabilities.sum()
        current_node = torch.multinomial(probabilities, 1).item()
        walk.append(current_node)

    return walk


def train_by_randomwalk(total_loss:float,edge_list:list,model:SkipGram,
                        loss_function=nn.CrossEntropyLoss(),label=0,walk_length=5):
    #model = SkipGram(num_nodes, emb_size, num_features)
    optimizer = optim.SGD(model.parameters(), lr=learning_rate)

    for edge in edge_list:
        src_node, dest_node, edge_weight, edge_features = edge           
        src_node_vec = Variable(torch.LongTensor([src_node]))
        src_node_features = node_features[src_node]
        context_nodes = random_walk(edge_list, src_node, walk_length)
        for context_node in context_nodes:
            context_node_vec = Variable(torch.LongTensor([context_node]))
            context_node_features = node_features[context_node]

            optimizer.zero_grad()
            scores = model.forward(src_node_vec,context_node_vec,src_node_features,context_node_features)
            
            if label > 0:
                src = torch.ones(1, 1)  # Positive label
                loss = loss_function.forward(scores.view(1,-1), src)
            elif label < 0:
                src = torch.zeros(1, 1)  # Negative label
                loss = loss_function.forward(scores.view(1,-1), src)
            else:
                loss = loss_function.forward(scores, src_node_vec)

            loss.backward()
            optimizer.step()
            total_loss += loss.item()

    return total_loss


def find_paths(edge_list, start_node, end_node, length, current_path=[]):
    current_path = current_path + [start_node]
    if start_node == end_node and len(current_path) == length:
        return [current_path]
    if len(current_path) > length:
        return []
    paths = []
    for edge in edge_list:
        src, dest, _, _ = edge
        if src == start_node:
            new_paths = find_paths(edge_list, dest, end_node, length, current_path)
            for new_path in new_paths:
                paths.append(new_path)
    return paths


def all_paths(edge_list:list,nodes:list,max_path_len=3):
    # Generate all paths of length 5
    all_paths = []
    for i in nodes:
        for j in nodes:
            if i != j:
                all_paths.append([i] + find_paths(edge_list, i, j, max_path_len))

    # Flatten the list of paths
    all_paths = [item for sublist in all_paths for item in sublist]
    return all_paths


def train_by_fullpath(model:SkipGram,all_paths:list, total_loss:float,
                      loss_function=nn.CrossEntropyLoss(),label=0):
    
    optimizer = optim.SGD(model.parameters(), lr=learning_rate)
    for path in all_paths:
        src_node = path[0]
        src_node_vec = torch.LongTensor(src_node)
        src_node_features = node_features[src_node]
        for context_node in path:
            context_node_vec = torch.LongTensor([context_node])
            context_node_features = node_features[context_node]

            optimizer.zero_grad()
            scores = model.forward(src_node_vec, context_node_vec, src_node_features, context_node_features)
            
            if label > 0:
                src_node = torch.ones(1, 1)  # Positive label
            elif label < 0:
                src_node = torch.zeros(1, 1)  # Negative label
            
            loss = loss_function.forward(scores.view(1, -1), src_node_vec)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()

    return total_loss
    

if __name__ == "__main__":
    input_graph = 'golimumab network'
    dpd_graph = load_sample_resnet(input_graph)
    required_nodetypes = ['Disease','Protein','FunctionalClass','Complex','SmallMol']
    dpd_components = dpd_graph.components(required_nodetypes) # to make sure graph is fully connected
    if isinstance(dpd_components[0],ResnetGraph):        
        edge_list,node_features,idx2objs,node_stats,edge_stats = dpd_components[0].rn2tensor()
        #num_nodes = max(max(e[0] for e in edge_list), max(e[1] for e in edge_list)) + 1
        num_nodes = len(idx2objs)
        num_features = len(node_features[0])
        
        training_set = dpd_components[0].subgraph_by_relprops(['ClinicalTrial'])
        training_edge_list,_,_,_,_ = training_set.rn2tensor(node_stats,node_features,edge_stats,idx2objs)
        training_node_idx = {e[0] for e in training_edge_list}
        training_node_idx.update({e[1] for e in training_edge_list})
        training_paths = all_paths(training_edge_list,list(training_node_idx))
        indications = training_set._psobjs_with('Disease','ObjTypeName')

        drugs = set(dpd_components[0]._psobjs_with('SmallMol','ObjTypeName'))
        diseases = dpd_components[0]._psobjs_with('Disease','ObjTypeName')
        diseases = [d for d in diseases if d not in indications]
        # to exclude indications reported as toxicities in clinical trials
        negative_set = dpd_components[0].neighborhood(drugs,diseases,['Regulation'],['positive'])
        negative_edge_list,_,_,_,_ = negative_set.rn2tensor(node_stats,node_features,edge_stats,idx2objs)
        neg_node_idx = {e[0] for e in training_edge_list}
        neg_node_idx.update({e[1] for e in training_edge_list})
        negative_paths = all_paths(negative_edge_list,list(neg_node_idx))

        test_set = dpd_components[0].neighborhood(drugs,diseases,['Regulation'],['negative'])
        test_set = test_set.subgraph_by_relprops(['Regulation'])
        test_edge_list,_,_,_,_ = test_set.rn2tensor(node_stats,node_features,edge_stats,idx2objs)
        
        # Hyperparameters
        emb_size = 2
        learning_rate = 0.01
        epochs = 1000

        # Create SkipGram model
        model = SkipGram(num_nodes, emb_size, num_features)
        #loss_function = nn.CrossEntropyLoss()
        loss_function = nn.BCEWithLogitsLoss()
        optimizer = optim.SGD(model.parameters(), lr=learning_rate)

        negative_samples = list(negative_set.edges())

        # Training loop
        for epoch in range(epochs):
            total_loss = 0.0
            total_loss = train_by_randomwalk(total_loss,training_edge_list,model,loss_function,label=1,walk_length=5)
            total_loss = train_by_randomwalk(total_loss,negative_edge_list,model,loss_function,label=-1,walk_length=5)

            if epoch % 100 == 0:
                print(f'Epoch: {epoch}, Loss: {total_loss}')

        # Get the learned embeddings
        print('Node embeddings:')
        print (model.node_embeddings.weight.data.numpy())
        print('Node feature embeddings:')
        print(model.feature_embeddings.weight.data.numpy())

        # Perform link prediction on a new pair of nodes
        with open(f'Link prediction in {input_graph}.txt','w',encoding='utf-8') as f:
            f.write('Drug\tDisease\tProbability\tRefCount\n')
            for drug_idx,disease_idx,_,_ in test_edge_list:
                drug = idx2objs[drug_idx]
                drug_name = drug.name()
                drug_urn = drug.urn()[:-1]
                drug_uid = PSObject.urn2uid(drug_urn)
                disease = idx2objs[disease_idx]
                disease_name = disease.name()
                disease_urn = disease.urn()[:-1]
                disease_uid = PSObject.urn2uid(disease_urn)
                rels = [rel for ruid, tuid, rel in dpd_components[0].edges.data('relation') 
                if ruid == drug_uid and tuid == disease_uid]
                refcount = 0
                for r in rels:
                    refcount += r.get_reference_count()

                drug_vec = torch.LongTensor([drug_idx])
                disease_vec = torch.LongTensor([disease_idx])
                drug_feature = torch.LongTensor(node_features[drug_idx])
                disease_feature = torch.LongTensor(node_features[disease_idx])

                scores = model.forward(drug_vec, disease_vec, drug_feature, disease_feature)
                probability = torch.sigmoid(scores).item()

                f.write(f'{drug_name}\t{disease_name}\t{round(probability,6)}\t{refcount}\n')