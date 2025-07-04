fn sckmer() {
    let mut targeted_taxids: Vec<&[u8]>;
    if lca.is_some() {
        // Parse set of desired taxonomic ranks
        let mut reports = kreports.iter().collect::<Vec<_>>();
        if let Some(ranks) = lca {
            let ranks_sets = ranks
                .iter()
                .map(|x| x.as_bytes())
                .collect::<HashSet<&[u8]>>();
            reports = reports
                .into_iter()
                .filter(|kr| ranks_sets.contains(kr.rank.as_slice()))
                .collect();
        }
        targeted_taxids = reports.into_iter().map(|kr| kr.taxid.as_slice()).collect();
        targeted_taxids = targeted_taxids
            .into_iter()
            .filter_map(|taxid| taxid_to_descendants.get(taxid))
            .flatten()
            .copied()
            .collect()
    } else {
        targeted_taxids = kreports.iter().map(|kr| kr.taxid.as_slice()).collect();
    }
}
