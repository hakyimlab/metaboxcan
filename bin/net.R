#! /usr/bin/env Rscript
set.seed(6767)
# functions
acat = function(pvec)
{
  TT = sum( tan( (0.5 - pvec) *pi ) )
  .5 - atan(TT / length(pvec)) / pi
}

get_community <- function(ss,m.network){
  # generate a graph
  g <- graph_from_data_frame(m.network, directed = FALSE)

  #Get all second degree nodes
  sel.nodes <- ego(g,order = 2,nodes = as.vector(ss$gene))

  # get the subgraph
  ff = data.frame(id = unlist(sel.nodes))
  g2 <- induced_subgraph(g, ff$id)

  # Detect community using second degree nodes
  lc <- cluster_louvain(g)

  return(lc)
}

plt_edges <- function(ss,m.network) {
  # ss: list containing significant genes
  # m.network: the interaction network
  g <- graph_from_data_frame(m.network, directed = FALSE)

  # filter the nodes to plot differentl for mets and genes
  mets <- ss %>% filter(!startsWith(gene,"E")) %>% pull(gene)
  genes <- ss %>% filter(startsWith(gene,"E")) %>% pull(gene)

  # select neighbours
  sel.nodes.m <- ego(g,order = 2,nodes = mets)
  sel.nodes.g <- ego(g,order = 1,nodes = genes)

  # get the subgraph
  ff = data.frame(id = unlist(sel.nodes.m))
  ff <- rbind(ff, data.frame(id = unlist(sel.nodes.g)))
  g3 <- induced_subgraph(g, ff$id)

  # get the edges df
  n.edges <- as.data.frame(get.edgelist(g3)) %>%
    setnames(.,colnames(.), c("from","to")) %>%
    left_join(m.network, by = c("from" = "source", "to" = "target")) %>%
    dplyr::rename(width = score)

  return(n.edges)

}

plt_nodes <- function(lc,scores,gene_info,met_info,multixcan){
  # get the nodes df
  lc_df <- data.frame(as.list(membership(lc))) %>% t() %>%
    data.frame() %>% setnames(., names(.), "group") %>%
    rownames_to_column(var = "label")

  n.nodes <- scores %>% dplyr::rename(label = gene) %>%
    mutate(value = zscore ^ 2) %>%
    inner_join(lc_df, by = "label") %>%
    left_join(gene_info %>% dplyr::select(ensembl_gene_id,preferred_name,annotation),
              by = c("label" = "ensembl_gene_id")) %>%
    left_join(met_info %>%
                dplyr::select(Compound_name,COMP_ID,HMDB,super_class,sub_class,
                              SUPER_PATHWAY,SUB_PATHWAY,description,Metabolite.ID) %>%
                mutate(COMP_ID = paste0("M", COMP_ID)),
              by = c("label" = "Metabolite.ID")) %>%
    left_join(multixcan %>% dplyr::select(gene,gene_name) %>% 
                dplyr::mutate(gene = gsub("\\..*","",gene)),
              by = c("label" = "gene")) %>% 
    dplyr::mutate(preferred_name = ifelse((is.na(preferred_name) & !is.na(gene_name)), gene_name, preferred_name)) %>% 
    dplyr::mutate(name = ifelse(is.na(Compound_name),preferred_name,Compound_name)) %>%
    dplyr::mutate(name = ifelse(is.na(name),label,name)) %>%
    dplyr::mutate(name = ifelse(startsWith(name,"ENSG0"),"",name)) %>%
    dplyr::mutate(annotation = ifelse(is.na(annotation),"Gene type: LncRNA",annotation)) %>%
    dplyr::mutate(desc = ifelse(!is.na(HMDB) & HMDB != "", glue::glue("HMDB: {HMDB}"),"")) %>%
    dplyr::mutate(desc = ifelse(!is.na(SUPER_PATHWAY) & SUPER_PATHWAY != "",
                         ifelse(desc == "",glue::glue("SUPER PATHWAY: {SUPER_PATHWAY}"),
                         glue::glue("{desc}<br>SUPER PATHWAY: {SUPER_PATHWAY}")),desc)) %>%
    dplyr::mutate(desc = ifelse((is.na(SUPER_PATHWAY) | SUPER_PATHWAY == "") & (!is.na(super_class) & super_class != ""),
                         ifelse(desc == "",glue::glue("SUPER CLASS: {super_class}"),
                         glue::glue("{desc}<br>SUPER CLASS: {super_class}")),desc)) %>%
    dplyr::mutate(desc = ifelse(!is.na(SUB_PATHWAY) & SUB_PATHWAY != "",
                         ifelse(desc == "",glue::glue("SUB PATHWAY: {SUB_PATHWAY}"),
                         glue::glue("{desc}<br>SUB PATHWAY: {SUB_PATHWAY}")),desc)) %>%
    dplyr::mutate(desc = ifelse((is.na(SUB_PATHWAY) | SUB_PATHWAY == "") &(!is.na(sub_class) & sub_class != ""),
                         ifelse(desc == "",glue::glue("SUB CLASS: {sub_class}"),
                         glue::glue("{desc}<br>SUB CLASS: {sub_class}")),desc)) %>%
    dplyr::mutate(desc = ifelse(desc == "" & !str_starts(label, "ENSG0"),
                         glue::glue("name: {name}"),desc)) %>%
    dplyr::mutate(desc = ifelse(!is.na(description) & description != "",
                         ifelse(desc == "",description,
                         glue::glue("{desc}<br>{description}")),desc)) %>%
    dplyr::mutate(desc = ifelse(desc == "",annotation,desc)) %>%
    dplyr::mutate(desc = glue::glue("ID: {label}<br>{desc}")) %>%
    dplyr::mutate(shape = ifelse(!startsWith(label,"ENSG0"), "square","dot")) %>%
    mutate(shadow = ifelse(b_sig == "yes", TRUE,FALSE)) %>%
    mutate(p2 = format(pvalue,scientific = T, digits = 2)) %>%
    mutate(desc = glue::glue("<p style='color: black; font-size: 10px;
                             font-family: Tahoma, sans-serif;
                             '><b>p-value: {p2}, BF sig: {b_sig}</b><br>{desc}</p>")) %>%
    #filter((label %in% n.edges$from) | (label %in% n.edges$to)) %>%
    dplyr::rename(id=label,label=name,title = desc) %>%
    dplyr::select(id,label,value,group,shape,title,shadow,pvalue,zscore)

  return(n.nodes)
}

get_neighbors <- function(m.network,ss,deg=1) {
  # find all the 1st degree neighbors for the significant nodes
  #TODO: Change gene id to gene names
  # generate a graph
  g <- graph_from_data_frame(m.network, directed = FALSE)

  #Get all second degree nodes
  sel.nodes <- ego(g,order = deg,nodes = as.vector(ss$gene))

  # get nodes 1st deg only
  vv = plyr::rbind.fill(lapply(sel.nodes, function(x)as.data.frame(t(names(x)))))

  df = vv %>% dplyr::rename(node = V1) %>%
    pivot_longer( cols = -node, names_to = "count", values_to = "members") %>%
    dplyr::select(-count) %>%
    dplyr::filter(!is.na(members)) %>%
    group_by(node) %>%
    summarize(total_n = n(),
      neighbors = paste(sort(unique(members)),collapse=", "))

  df <- ss %>% dplyr::select(-b_sig) %>%
    inner_join(df, by = c("gene" = "node")) %>% dplyr::rename(sig_node = gene)

  return(df)
}

get_clusters <- function(nodes){

  cluster_summary <- nodes %>% group_by(group) %>%
    arrange(pvalue) %>%
    transmute(mean_zscore = mean(zscore),
              min_zscore = min(zscore),max_zscore = max(zscore),
              acat_pvalue = acat(pvalue),
              no_of_members = n(),
              neighbors = paste(label[1:10],collapse=", ")) %>%
    distinct() %>% mutate(group = paste0("cluster_",group)) %>%
    dplyr::mutate(neighbors = str_remove_all(neighbors, ", NA")) %>%
    mutate_at(vars(mean_zscore, min_zscore, max_zscore), ~round(., 2)) %>%
    arrange(acat_pvalue) %>%
    mutate(acat_pvalue = format(acat_pvalue, scientific= T, digits=2))
  # return
  return(cluster_summary)
}


top_selector <- function(scores,network) {
  # purpose is to select manageable number of nodes to visualize
  scores <- scores %>%
    filter((gene %in% network$source) |(gene %in% network$target))
  ## selecting top 100 nodes
  if(sum(scores$b_sig == "yes", na.rm = T) < 100){
    pp <- 95 - sum(scores$b_sig == "yes", na.rm = T)
    # slice up
    ss <- scores %>% filter(b_sig == "yes") %>%
      bind_rows(scores %>% filter(b_sig == "no") %>%
                  slice_min(pvalue, n = pp))
  } else {
    ss <- scores %>% filter(b_sig == "yes")
  }
  return(ss)
}

