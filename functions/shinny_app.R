get_click_data = function(x) {
  with(x, {
    tryCatch({
      
      case_id = strsplit(domain$discrete_limits$x[[round(x)]], "/")[[1]]
      
      list(
        panel_x = panelvar1,
        panel_y = panelvar2,
        case = case_id[1],
        tissue = case_id[2],
        peak = domain$discrete_limits$y[[round(y)]]
      )
      
    }, error = function(e)
      return(NULL))
  })
}



get_stat_text = function(idx) {
  tryCatch({
    # insert column index
    idx_cn = paste0(idx$case, "-normal_vs_", tolower(idx$tissue))
    mt = match(tolower(idx_cn), tolower(colnames(final_rec_count$p)))
    p_value_edger = final_rec_count$p[idx$peak, mt]
    effect_size = final_rec_count$fc[idx$peak, mt]
    p_value_deseq = NA
    p_value_deseq_sc = NA
    
    paste0(
      "Peak: ", idx$peak, "\n",
      "Case: ", idx$case, "\n",
      "edgeR model: \n",
      "  - p-value: ", format.pval(p_value_edger), "\n",
      "  - effect size: ", effect_size, "\n",
      "- p-value of DESeq2 model: ", p_value_deseq, "\n",
      "- p-value of DESeq2 model (subclonal): ", p_value_deseq_sc, "\n"
    )
  }, error = function(e)
    return(""))
}

