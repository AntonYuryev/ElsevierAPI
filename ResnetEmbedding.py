# coding=utf-8
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
from tf_geometric.utils import tf_utils
import tf_geometric as tfg
from tf_geometric.data.graph import Graph as tfgraph
import tensorflow as tf
from tf_geometric.utils.graph_utils import edge_train_test_split, negative_sampling
from tqdm import tqdm
from contextlib import redirect_stdout
from ElsevierAPI.ResnetAPI.ResnetAPIcache import APIcache, PSObject, PSRelation
from ElsevierAPI.pandas.panda_tricks import df,pd,np
from ElsevierAPI.ResnetAPI.ResnetGraph import ResnetGraph, REGULATORS, TARGETS


def load_resnet():
    prot2disease = APIcache(cache_name='prot_reg_disease',connect2server=False)
    drug2prot = APIcache(cache_name='drug2target',connect2server=False)
    pppi = APIcache(cache_name='PPPI network',connect2server=False)    
    resnet_graph = prot2disease.network.compose(drug2prot.network)
    resnet_graph = resnet_graph.compose(pppi.network)
    print(f'Loaded resnet graph with {len(resnet_graph)} nodes and {resnet_graph.number_of_edges()} relations')
    return resnet_graph


def load_sample_resnet(cache_name:str):
    sample = APIcache(cache_name=cache_name,connect2server=False)
    return sample.network


def set_negative_edges(resnet:ResnetGraph,max_num=100):
    my_graph = resnet.copy()
    drugs = my_graph._psobjs_with('SmallMol','ObjTypeName')
    diseases = my_graph._psobjs_with('Disease','ObjTypeName')
    edge_counter = 0
    for drug in drugs:
        if edge_counter >= max_num: break
        for disease in diseases:
            drug_uid = drug.uid()
            disease_uid = disease.uid()
            if not my_graph.has_edge(drug_uid,disease_uid):
                new_rel = PSRelation({'ObjTypeName':['EdgeNotExist']})
                new_rel.Nodes[REGULATORS] = [(drug_uid, '0', 'unknown')]
                new_rel.Nodes[TARGETS] = [(disease_uid, '1', 'unknown')]
                edge_counter += 1
                if edge_counter < max_num:
                    my_graph.add_rel(new_rel)
                else:
                    break
    return my_graph


def resnet2tfg(resnet:ResnetGraph):
    node_pd = pd.DataFrame(columns=['URN','Type'])
    train_edge_rows = list()
    test_neg_edge_rows = list()
    test_edge_rows  = list()
    
    my_graph = set_negative_edges(resnet)

    def node_idx(n:PSObject,in_df:pd.DataFrame):
        pdcopy = in_df.copy()
        urn = n.urn()
        activated_urn = urn + 'a'
        repressed_urn = urn + 'i'

        acitvated_idx = pdcopy[pdcopy['URN'] == activated_urn].index.tolist()
        if not acitvated_idx:
            pdcopy = pd.concat([pdcopy,pd.DataFrame({'URN':[activated_urn], 'Type':[n.objtype()]})], ignore_index=True)
        
        repressed_idx = pdcopy[pdcopy['URN'] == repressed_urn].index.tolist()
        if not repressed_idx:
           pdcopy = pd.concat([pdcopy,pd.DataFrame({'URN':[repressed_urn], 'Type':[n.objtype()]})], ignore_index=True)

        return pdcopy[pdcopy['URN'] == activated_urn].index.tolist(), pdcopy[pdcopy['URN'] == repressed_urn].index.tolist(), pdcopy


    edge_rows = list()
    rel_idx = 0
    for r,t,urn,e in my_graph.edges.data(keys=True):
        psrel = e['relation']
        if isinstance(psrel,PSRelation):
            regulator = my_graph._get_node(r)
            target = my_graph._get_node(t)

            active_regulator_idx,repressed_regulator_idx,node_pd = node_idx(regulator,node_pd)
            active_target_idx,repressed_target_idx,node_pd = node_idx(target,node_pd)
            reltype = psrel.objtype()
            if reltype == 'ClinicalTrial':
                effect_sign = -1 # hacking ClinicalTrial that do not have Effect sign in Resnet
            else:
                effect_sign = psrel.effect_sign()

            relweight = psrel.get_reference_count()
            
            if effect_sign > 0:
                my_rt_tup = [(active_regulator_idx[0],active_target_idx[0],relweight,reltype),
                                (repressed_regulator_idx[0],repressed_target_idx[0],relweight,reltype)]
            elif effect_sign < 0:
                my_rt_tup = [(active_regulator_idx[0],repressed_target_idx[0],relweight,reltype),
                                (repressed_regulator_idx[0],active_target_idx[0],relweight,reltype)]
            else:
                my_rt_tup = [(active_regulator_idx[0],active_target_idx[0],relweight,reltype),
                                (repressed_regulator_idx[0],repressed_target_idx[0],relweight,reltype),
                                (repressed_regulator_idx[0],active_target_idx[0],relweight,reltype),
                                (active_regulator_idx[0],repressed_target_idx[0],relweight,reltype)]

            if reltype != 'EdgeNotExist':
                edge_rows += my_rt_tup

            if regulator.objtype() == 'SmallMol':
                if target.objtype() == 'Disease':
                    if reltype == 'ClinicalTrial':
                        train_edge_rows += my_rt_tup
                    elif reltype == 'Regulation':
                        test_edge_rows += my_rt_tup
                    else:
                        assert(reltype == 'EdgeNotExist')
                        test_neg_edge_rows += my_rt_tup

            rel_idx += 1
    
    edge_df = df.from_rows(edge_rows,header=['Regulator_idx','Target_idx','RefCount','RelType'])
    train_edge_df =  df.from_rows(train_edge_rows,header=['Regulator_idx','Target_idx','RefCount','RelType'])
    test_edge_df =  df.from_rows(test_edge_rows,header=['Regulator_idx','Target_idx','RefCount','RelType'])
    test_neg_edge_df =  df.from_rows(test_neg_edge_rows,header=['Regulator_idx','Target_idx','RefCount','RelType'])
    
    one_hot_encoded_data = pd.get_dummies(node_pd, columns = ['Type'])
    one_hot_encoded_data = one_hot_encoded_data.drop(columns=['URN'])
    x = np.ascontiguousarray(one_hot_encoded_data.to_numpy(float))

    edge_index = edge_df[['Regulator_idx','Target_idx']].transpose().to_numpy(int)
    train_edge_index = train_edge_df[['Regulator_idx','Target_idx']].transpose().to_numpy(int)
    test_edge_index = test_edge_df[['Regulator_idx','Target_idx']].transpose().to_numpy(int)
    test_neg_edge_index = test_neg_edge_df[['Regulator_idx','Target_idx']].transpose().to_numpy(int)

    edge_df = edge_df.l2norm({'RefCount':'weight'})
    edge_weight = edge_df['weight'].transpose().to_numpy(float)
    train_edge_df = train_edge_df.l2norm({'RefCount':'weight'})
    train_edge_weight = train_edge_df['weight'].transpose().to_numpy(float)

    edge_type = edge_df['RelType'].transpose() # to do: encode edge types here

    #param x: Tensor/NDArray, shape: [num_nodes, num_features], node features
    #:param edge_index: Tensor/NDArray, shape: [2, num_edges], edge information.
    #Each column of edge_index (u, v) represents an directed edge from u to v.
    #:param edge_weight: Tensor/NDArray/None, shape: [num_edges]
    graph = tfgraph(x,edge_index,edge_weight=edge_weight)
    train_graph = tfgraph(x,train_edge_index,edge_weight=train_edge_weight)
    print (f'Loaded graph containing {graph.num_nodes} nodes with {graph.num_features} features,\
and {graph.num_edges} edges')
    
    print (f'Training graph contains {train_graph.num_nodes} nodes with {train_graph.num_features} features,\
and {train_graph.num_edges} ClinicalTrial edges from Resnet')
    
    print (f'Dataset has {len(set(test_edge_index[1]))/2} diseases for postive test. These diseases are linked to drug by Regulation in Resnet')
    print (f'Dataset has {len(set(test_neg_edge_index[1]))/2} diseases for negative test. These diseases are not liked to drug in Resnet.')
    
    return graph, train_graph, test_edge_index, test_neg_edge_index


@tf_utils.function
def encode(graph:tfgraph, training=False):
    h = gcn0([graph.x, graph.edge_index, graph.edge_weight], cache=graph.cache) #shape=(2708, 32)
    #print(h)
    h = dropout(h, training=training)
    h = gcn1([h, graph.edge_index, graph.edge_weight], cache=graph.cache) #shape=(2708, 16)
    return h


@tf_utils.function
def predict_edge(embedded, edge_index):
    row, col = edge_index[0], edge_index[1]
    embedded_row = tf.gather(embedded, row)
    embedded_col = tf.gather(embedded, col)

    # dot product
    logits = tf.reduce_sum(embedded_row * embedded_col, axis=-1)
    return logits


@tf_utils.function
def compute_loss(pos_edge_logits, neg_edge_logits):
    pos_losses = tf.nn.sigmoid_cross_entropy_with_logits(
        logits=pos_edge_logits,
        labels=tf.ones_like(pos_edge_logits)
    )

    neg_losses = tf.nn.sigmoid_cross_entropy_with_logits(
        logits=neg_edge_logits,
        labels=tf.zeros_like(neg_edge_logits)
    )

    return tf.reduce_mean(pos_losses) + tf.reduce_mean(neg_losses)


def evaluate(_test_edge_index,_test_neg_edge_index):
    embedded = encode(train_graph)
    pos_edge_logits = predict_edge(embedded, _test_edge_index)
    neg_edge_logits = predict_edge(embedded, _test_neg_edge_index)

    pos_edge_scores = tf.nn.sigmoid(pos_edge_logits)
    print('pos_edge_scores in test_edge_index:')
    print(pos_edge_scores)
    neg_edge_scores = tf.nn.sigmoid(neg_edge_logits)
    print('neg_edge_scores in test_neg_edge_index:')
    print(neg_edge_logits)

    y_true = tf.concat([tf.ones_like(pos_edge_scores), tf.zeros_like(neg_edge_scores)], axis=0)
    y_pred = tf.concat([pos_edge_scores, neg_edge_scores], axis=0)

    auc_m = tf.keras.metrics.AUC()
    auc_m.update_state(y_true, y_pred)

    return auc_m.result().numpy()


if __name__ == "__main__":
    # dpd_graph = load_resnet() # dpd = disease-protein-drug
    dpd_graph = load_sample_resnet('golimumab network')
    required_nodetypes = ['Disease','Protein','FunctionalClass','Complex','SmallMol']
    dpd_components = dpd_graph.components(required_nodetypes) # to make sure we have fully connected graph
    graph, train_graph, test_edge_index, test_neg_edge_index = resnet2tfg(dpd_components[0])

    with open('LinkPredict.log', 'w') as fout:
        with redirect_stdout(fout):
            embedding_size = 16
            drop_rate = 0.2

            gcn0 = tfg.layers.GCN(32, activation=tf.nn.relu)
            gcn1 = tfg.layers.GCN(embedding_size)
            dropout = tf.keras.layers.Dropout(drop_rate)

            gcn0.build_cache_for_graph(graph)
            gcn0.build_cache_for_graph(train_graph)

            optimizer = tf.keras.optimizers.Adam(learning_rate=1e-2)

            for step in tqdm(range(1000)):
                with tf.GradientTape() as tape:
                    embedded = encode(train_graph, training=True)

                    # negative sampling for training
                    train_neg_edge_index = negative_sampling(
                        train_graph.num_edges,
                        graph.num_nodes,
                        edge_index=None#train_graph.edge_index
                    )

                    pos_edge_logits = predict_edge(embedded, train_graph.edge_index)
                    neg_edge_logits = predict_edge(embedded, train_neg_edge_index)

                    loss = compute_loss(pos_edge_logits, neg_edge_logits)

                vars = tape.watched_variables()
                grads = tape.gradient(loss, vars)
                optimizer.apply_gradients(zip(grads, vars))

                if step % 20 == 0:
                    auc_score = evaluate(test_edge_index, test_neg_edge_index)
                    print("step = {}\tloss = {}\tauc_score = {}".format(step, loss, auc_score))
