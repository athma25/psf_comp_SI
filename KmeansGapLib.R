library(tidyverse)
library(FactoClass)

## Calculate within-cluster sum of squares (wcss) for a given community and number of clusters
kmeans_withinss = 
  function(traits, abundances, numclusters, null, nstart){
    
    set.seed(null)
    
    if (null != 0){
			traits = runif(length(traits), min = min(traits), max = max(traits))
			abundances = sample(abundances)
		}
 
		kmw =
      kmeansW(
        x = traits,
        centers = numclusters,
        nstart = nstart,
        weight = abundances
      )
   
    wcss =
      kmw %>%
      pluck('withinss') %>%
      mean()
	
		centers =
			kmw %>%
			pluck('centers')

    return(lst(kmean_sum=tibble(clusters = numclusters, null = null, wcss = wcss),centers=tibble(centers)))
  }

## Catalog wcss for all null randomizations of the community (including the original), up to maxnumclusters 
kmeans_prep = 
  function(data, trait_index, abundance_index, maxnumclusters, numnulls, nstart){
    
    traits = data[, trait_index]
    abundances = data[, abundance_index]
    
    parms = 
      expand_grid(
        numclusters = 2:maxnumclusters,
        null = 0:numnulls
      )
    
    foo = function(numclusters, null) 
      kmeans_withinss(traits, abundances, numclusters, null, nstart)$kmean_sum
    
    res = 
      parms %>%
      pmap_dfr(foo)
    
    return(res)
  }


## Kmeansgap analysis on run# from batch id
kmeans_gap =
	function(id, run, rep){
	
		fOut=sprintf("output/%s/Kmeans_summary.txt",id);
		fOut_ls=sprintf("output/%s/ls_summary.txt",id);
	
		## Fixed parameters
		number_of_seeds = 100
		number_of_nulls = 100
		max_number_of_clusters = 25
		trait_index = 1
		abundance_index = 7

		for(i in 0:(rep-1)){
			filename=sprintf("output/%s/run%03d_%03d_final.txt",id,run,i);
			fPlot=sprintf("output/%s/run%03d_%03d_kmeans.png",id,run,i);

			## Read datafile
			data = read.table(filename, sep = ',')
			data = data[-1,]			### Remove first row which contains run summary

			## Generate tibble by appling kmeans_prep on data
			dat = 
			kmeans_prep(
				data = data, 
				trait_index = trait_index, 
				abundance_index = abundance_index, 
				maxnumclusters = max_number_of_clusters, 
				numnulls= number_of_nulls,
				nstart = number_of_seeds
			)

			## Summary statistics (mean and sd) for each number of clusters across nulls
			stats = 
				dat %>% 
				group_by(clusters) %>% 
				summarize(mean_wcss = mean(wcss), sd_wcss = sd(wcss))

			## Create tibble showing the gap in log(wcss) for each number of clusters and each community (original and nulls)
			gap =
				stats %>%
				inner_join(
					dat,
					by = 'clusters'
				) %>%
				mutate(gap = log(mean_wcss) - log(wcss))

			## Find 95% quantile of the maxgap statistic across nulls 
			## (if peak is above the red line we have statisitcal significance)
			quants =
				gap %>%
				group_by(null) %>%
				slice_max(gap) %>%
				ungroup() %>%
				summarize(qtl = quantile(gap, .95))


			## P-value
			nullgap = 
				gap %>% 
				group_by(null) %>% 
				slice_max(gap) %>%
				pull(gap)

			thegap = 
				gap %>% 
				filter(null == 0) %>% 
				slice_max(gap) %>% 
				pull(gap)

			pval = mean(thegap < nullgap)


			## Log result
			result = 
				gap %>%
				filter(null == 0) %>%
				slice_max(gap) %>%
				bind_cols(quants) %>%
				bind_cols(p.value = pval)
		
			write_delim(bind_cols(the_run=run,the_rep=i,result),fOut,col_names=FALSE,append=TRUE);

			png(file=fPlot);
			## Plot gap curve
			plot_gap = 
				gap %>%
				filter(null == 0) %>%
				ggplot(aes(clusters, gap)) +
				geom_line() +
				geom_point() +
				geom_hline(
					aes(yintercept = qtl), 
					data = quants, 
					color = 'red'
				)

			## Plot stems
			plot_stems = 
				data %>%
				as_tibble() %>%
				ggplot() +
				geom_segment(aes(x = V1, xend = V1, y = rep(0, nrow(data)), yend = V7)) +
				geom_point(aes(V1, V7)) +
				labs(x = 'trait', y = 'abundance')

			gridExtra::grid.arrange(plot_stems, plot_gap)
			dev.off();

			## Limiting similarity
			centers=kmeans_withinss(data[, trait_index],data[, abundance_index],result$clusters,0,number_of_seeds)$centers
			if (nrow(centers)>2){
				centers=sort(centers$centers)
				var_centers=var(head(c(centers[-1],NaN)-centers,-1))

				rand_centers=matrix(runif(number_of_nulls*result$clusters,min=min(data[,trait_index]),max=max(data[,trait_index])),nrow=result$clusters);
				rand_centers=apply(rand_centers,2,sort);
				null_var=apply(head(rbind(rand_centers[-1,],NaN)-rand_centers,-1),2,var);
	
				ls_quants=quantile(null_var,0.05)
	
				pval=mean(var_centers > null_var)
	
				ls_result=bind_cols(the_run=run,the_rep=i,the_var=var_centers,ls_qtl=ls_quants,p.value=pval)
			} else {
				ls_result=bind_cols(the_run=run,the_rep=i,the_var=999,ls_quants=999,p.value=1)
				print("Method detected too few or too many clusters");
			}

			write_delim(ls_result,fOut_ls,col_names=FALSE,append=TRUE);

		print(sprintf("Completed run%03d_%03d",run,i));
		}
	}	
