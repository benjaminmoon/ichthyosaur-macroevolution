dist_bin <- lapply(seq_along(1:2), function (scheme) {
              lapply(seq_along(1:4), function (dist) {
                list(ddat = dist_data[[dist]],
                     binning = bin_data[[scheme]]$bin_data,
                     dist = plot_names$dist[dist],
                     bins = names(plot_names$bins[scheme]),
                     corr = "",
                     metric = "pd")
              })
            }) %>%
              unlist(recursive = FALSE)

pd_rare <- lapply(seq_along(dist_bin), function (run) {
             lapply(seq_along(dist_bin[run]$binning), function (bin) {
               # get number of rarefactions to do
               nrar <- seq(2, length(dist_bin[run]$binning[bin]))
           
               # sample taxa
               lapply(nrar, function (n_samp) {
                 replicate(n = 500,
                           sample(bin_data[run$dist]$bin_data[bin],
                                  size = n_samp,
                                  replace = FALSE))
               })
             })
           })
