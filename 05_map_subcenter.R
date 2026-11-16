#!/usr/bin/env Rscript
# Universal version of the Paris FAR/RCOI/C-I-O layout, with the same analytical rules.
suppressPackageStartupMessages({library(sf);library(dplyr);library(readr);library(ggplot2);library(patchwork)})
a<-grep('^--file=',commandArgs(),value=TRUE)
root<-normalizePath(file.path(dirname(sub('^--file=','',a[1])),'../..'));setwd(root)
input_dir<-'outputs/gt5_subcentres_far_rcoi_poi_details_redorange'
args<-commandArgs(trailingOnly=TRUE)
stopifnot(all(grepl('^--(cities|output|above-mean|minimal-text)=.+$',args)),!anyDuplicated(sub('=.*','',args)))
arg<-function(key,default){z<-args[grepl(paste0('^--',key,'='),args)];if(length(z))sub('^[^=]*=','',z)else default}
out<-arg('output',input_dir)
above_mean<-arg('above-mean','false');stopifnot(above_mean%in%c('true','false'));above_mean<-above_mean=='true'
minimal_text<-arg('minimal-text','false');stopifnot(minimal_text%in%c('true','false'));minimal_text<-minimal_text=='true'
if(above_mean && out==input_dir)stop('Use a separate --output directory for the above-mean filter')
for(d in c(out,file.path(out,'maps'),file.path(out,'gpkg_by_city')))dir.create(d,recursive=TRUE,showWarnings=FALSE)
eligible<-read_csv(file.path(input_dir,'eligible_cities.csv'),show_col_types=FALSE)
fit<-read_csv('outputs/far_power_law_41_footprint_cut/far_power_law_fit_41_cities.csv',show_col_types=FALSE)
cities<-fit$city[fit$city%in%eligible$city]
requested<-arg('cities','');if(nzchar(requested)){z<-strsplit(requested,',',fixed=TRUE)[[1]];stopifnot(all(z%in%cities),!anyDuplicated(z));cities<-z}
poi_manifest<-read_csv(file.path(input_dir,'poi_source_manifest.csv'),show_col_types=FALSE)
road_manifest<-read_csv(file.path(input_dir,'road_source_manifest.csv'),show_col_types=FALSE)
lookup<-read_csv('outputs/subcentre_standardized_poi_rcoi_32_cities/rcoi_16_type_colours.csv',show_col_types=FALSE)
lookup$display_type<-ifelse(lookup$outcome=='Mixed','RCOI',ifelse(lookup$outcome=='Non-subcentre','LLLL (all low)',lookup$outcome))
original_refs<-read_csv('outputs/subcentre_standardized_poi_rcoi_32_cities/city_standardization_reference.csv',show_col_types=FALSE)
source<-'outputs/subcentres_central_far_30_percent_41_cities/All_cities_central_FAR_30_percent.gpkg'
all_patches<-st_read(source,layer='threshold_patches',quiet=TRUE)
all_cells<-st_read(source,layer='selected_cells',quiet=TRUE)
centres<-st_read(source,layer='urban_centres',quiet=TRUE);bounds<-st_read(source,layer='insee_boundaries',quiet=TRUE)
far_source<-'outputs/all_city_cells_zipf_boundaries_500m_buffer_footprint_cut_reused_buffer_area/all_city_cells_zipf_boundaries_500m_buffer_footprint_cut_reused_buffer_area.gpkg'
far_levels<-c('<0.4','0.4–<0.8','0.8–<1','1–2','>2–4','>4','No valid FAR')
far_colours<-setNames(c('#FFF7EC','#FEE8C8','#FDD49E','#FDBB84','#EF6548','#B30000','#EEEEEE'),far_levels)
poi_names<-c('Commercial','Industrial','Office / employment')
poi_colours<-c('Commercial'='#E69F00','Industrial'='#7B3294','Office / employment'='#0072B2')
poi_shapes<-c('Commercial'=16,'Industrial'=17,'Office / employment'=0)
poi_sizes<-c('Commercial'=.85,'Industrial'=1.15,'Office / employment'=1.25)
metrics<-c(R='residence',C='commercial',O='employment',I='industrial')
slug<-function(x)gsub('[^A-Za-z0-9]+','_',stringi::stri_trans_general(x,'Latin-ASCII'))
legacy_slug<-function(x)gsub('[^A-Za-z0-9]+','_',iconv(x,to='ASCII//TRANSLIT'))
for(e in parse('scripts/subcentre_city_atlas/map_subcentres_types_poi_percentages.R'))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('north_direction')))eval(e)
theme_map<-theme_void(base_size=10)+theme(plot.title=element_text(face='bold',size=13,hjust=0),plot.subtitle=element_text(size=9,hjust=0),
 legend.position='bottom',legend.title=element_text(size=9,face='bold'),legend.text=element_text(size=8),plot.margin=margin(7,7,7,7))
adorn<-function(p,bb,scale_m,text_size=2.7){
 xx<-unname(bb[c('xmin','xmax')]);yy<-unname(bb[c('ymin','ymax')]);span<-diff(xx)
 sx<-xx[1]+span*.07;sy<-yy[1]+diff(yy)*.07;tick<-diff(yy)*.008
 ax<-xx[2]-span*.075;ay<-yy[2]-diff(yy)*.14;nd<-north_direction(ax,ay);perp<-c(-nd[2],nd[1]);al<-diff(yy)*.075
 tip<-c(ax,ay)+nd*al;tail<-c(ax,ay)-nd*al*.45;left<-c(ax,ay)+perp*al*.23;right<-c(ax,ay)-perp*al*.23
 p+annotate('segment',x=sx,xend=sx+scale_m,y=sy,yend=sy,linewidth=.6)+
  annotate('segment',x=c(sx,sx+scale_m),xend=c(sx,sx+scale_m),y=sy-tick,yend=sy+tick,linewidth=.3)+
  annotate('text',x=sx+scale_m/2,y=sy+tick*3,label=if(scale_m>=1000)paste(scale_m/1000,'km')else paste(scale_m,'m'),size=text_size)+
  annotate('polygon',x=c(tip[1],left[1],tail[1]),y=c(tip[2],left[2],tail[2]),fill='#111111',color='#111111',linewidth=.3)+
  annotate('polygon',x=c(tip[1],right[1],tail[1]),y=c(tip[2],right[2],tail[2]),fill='white',color='#111111',linewidth=.3)+
  annotate('text',x=ax+al*nd[1]*1.3,y=ay+al*nd[2]*1.3,label='N',fontface='bold',size=text_size+.3)
}
classify_patches<-function(city,patches,selected,cc){
 means<-vapply(metrics,function(k)mean(cc[[paste0(k,'_count')]]),numeric(1))
 sds<-vapply(metrics,function(k)sd(cc[[paste0(k,'_count')]]),numeric(1));stopifnot(all(sds>0))
 ref<-original_refs[original_refs$city==city,]
 if(nrow(ref)){ref<-ref[match(names(metrics),ref$initial),];stopifnot(max(abs(ref$city_mean_count-means))<1e-8,max(abs(ref$city_sd_count-sds))<1e-8)}
 at<-match(as.character(selected$id),as.character(cc$id))
 rows<-lapply(seq_len(nrow(patches)),function(j){
  hit<-at[selected$patch_id==patches$patch_id[j]];hit<-hit[!is.na(hit)]
  pm<-vapply(metrics,function(k)if(length(hit))mean(cc[[paste0(k,'_count')]][hit])else NA_real_,numeric(1))
  z<-(pm-means)/sds;code<-if(length(hit))paste0(ifelse(z>1e-10,'H','L'),collapse='')else NA_character_
  index<-match(code,lookup$code)
  row<-data.frame(patch_id=patches$patch_id[j],analysed_cells=length(hit),poi_cell_coverage_fraction=length(hit)/patches$cells[j],rcoi_code=code,
   rcoi_type=if(length(hit))lookup$display_type[index]else'POI unavailable',rcoi_colour=if(length(hit))lookup$type_colour[index]else'#999999')
  for(k in names(metrics)){row[[paste0(k,'_mean_count')]]<-pm[k];row[[paste0(k,'_mean_z')]]<-z[k]}
  row
 })
 list(patches=left_join(patches,bind_rows(rows),by='patch_id'),reference=data.frame(city=city,category=names(metrics),city_cells=nrow(cc),city_mean_count=means,city_sd_count=sds))
}
choose_windows<-function(patches,poi,viewport,ee){
 cp<-st_read(poi,query='SELECT osm_uid,geom FROM commercial',wkt_filter=st_as_text(viewport),quiet=TRUE)
 candidates<-patches[!patches$is_central_core & patches$analysed_cells>0 & lengths(st_intersects(patches,viewport))>0,]
 candidates$commercial_count<-lengths(st_intersects(candidates,cp))
 candidates<-candidates|>arrange(desc(commercial_count),desc(area_km2),patch_id)
 chosen<-integer();boxes<-list();overlap<-FALSE
 for(pass in 1:2){
  for(j in seq_len(nrow(candidates))){
   if(length(chosen)==3)break
   if(j%in%chosen)next
   xy<-unname(st_coordinates(suppressWarnings(st_point_on_surface(st_geometry(candidates[j,]))))[1,])
   bb<-c(xmin=xy[1]-1000,ymin=xy[2]-1000,xmax=xy[1]+1000,ymax=xy[2]+1000)
   if(bb['xmin']<ee['xmin']||bb['xmax']>ee['xmax']||bb['ymin']<ee['ymin']||bb['ymax']>ee['ymax'])next
   box<-st_as_sfc(st_bbox(bb,crs=3035))
   hits<-length(boxes)>0&&any(vapply(boxes,function(z)lengths(st_intersects(z,box))>0,logical(1)))
   if(pass==1&&hits)next
   chosen<-c(chosen,j);boxes[[length(boxes)+1]]<-box;overlap<-overlap||hits
  }
  if(length(chosen)==3)break
 }
 stopifnot(length(chosen)==3)
 target<-candidates[chosen,]
 list(windows=st_sf(detail_area=1:3,patch_id=target$patch_id,patch_label=target$patch_label,
   geometry=st_sfc(lapply(boxes,function(z)z[[1]]),crs=3035)),overlap=overlap)
}
get_context<-function(city,window){
 metadata<-file.path('outputs/subcentre_city_atlas_types_poi_percentages/building_cache',paste0(legacy_slug(city),'.rds'))
 if(!file.exists(metadata))metadata<-file.path('outputs/all_41_cities_functions_boundaries_all_buildings/building_cache',paste0(legacy_slug(city),'.rds'))
 bm<-readRDS(metadata)
 parts<-lapply(bm$signature,function(s)st_read(s$path,query=paste0('SELECT geom FROM ',s$layer),wkt_filter=st_as_text(st_geometry(window)),quiet=TRUE))
 buildings<-st_sf(geometry=do.call(c,lapply(parts,st_geometry)))
 sources<-road_manifest[road_manifest$city==city,];stopifnot(nrow(sources)>0)
 road_parts<-lapply(seq_len(nrow(sources)),function(j){
  p<-sources$road_path[j];layer<-sources$road_layer[j]
  field<-st_read(p,query=paste0("SELECT column_name FROM gpkg_geometry_columns WHERE table_name='",layer,"'"),quiet=TRUE)$column_name
  one<-st_read(p,query=paste0('SELECT ',field,' FROM ',layer,' LIMIT 1'),quiet=TRUE)
  w<-st_transform(window,st_crs(one))
  st_read(p,query=paste0('SELECT ',field,' FROM ',layer),wkt_filter=st_as_text(st_geometry(w)),quiet=TRUE)|>st_transform(3035)
 })
 roads<-st_sf(geometry=do.call(c,lapply(road_parts,st_geometry)))
 list(buildings=buildings,roads=roads)
}
make_city<-function(city){
 message('Mapping ',city)
 stem<-slug(city);gpkg<-file.path(out,'gpkg_by_city',paste0(stem,'_FAR_RCOI_CIO_min4.gpkg'))
 patches<-all_patches[all_patches$city==city & all_patches$cells>=4,]
 selected<-all_cells[all_cells$city==city & all_cells$patch_id%in%patches$patch_id,]
 secondary_count<-sum(!patches$is_central_core);stopifnot(secondary_count>5)
 pm<-poi_manifest[poi_manifest$city==city,];stopifnot(nrow(pm)==1)
 cc<-readRDS(pm$cell_cache)$cells
 classified<-classify_patches(city,patches,selected,cc);patches<-classified$patches
 excluded<-patches[FALSE,]
 if(above_mean){
  keep<-!is.na(patches$rcoi_code)&grepl('H',patches$rcoi_code)
  excluded<-patches[!keep,];patches<-patches[keep,]
  secondary_count<-sum(!patches$is_central_core)
 }
 centre<-centres[centres$city==city,];boundary<-bounds[bounds$city==city,]
 xy<-unname(st_coordinates(centre)[1,]);cx<-round(xy[1]/250)*250;cy<-round(xy[2]/250)*250
 ee<-c(xmin=cx-20000,ymin=cy-20000,xmax=cx+20000,ymax=cy+20000)
 viewport<-st_as_sfc(st_bbox(ee,crs=3035))
 chosen<-choose_windows(patches,pm$poi_path,viewport,ee);windows<-chosen$windows
 visible<-patches[lengths(st_intersects(patches,viewport))>0,]
 used<-c(lookup$display_type[lookup$display_type%in%visible$rcoi_type],if(any(visible$rcoi_type=='POI unavailable'))'POI unavailable')
 type_colours<-c(setNames(lookup$type_colour,lookup$display_type),'POI unavailable'='#999999')[used]
 far<-st_read(far_source,query=paste0("SELECT id,city,FAR,geom FROM all_city_cells_footprint_cut_reused_buffer_area WHERE city='",gsub("'","''",city),"'"),
   wkt_filter=st_as_text(viewport),quiet=TRUE)
 far$far_class<-factor(case_when(!is.finite(far$FAR)|far$FAR<0~far_levels[7],far$FAR<.4~far_levels[1],far$FAR<.8~far_levels[2],
  far$FAR<1~far_levels[3],far$FAR<=2~far_levels[4],far$FAR<=4~far_levels[5],TRUE~far_levels[6]),levels=far_levels)
 overview<-ggplot()+geom_sf(data=far,aes(fill=far_class),color=NA,show.legend=TRUE)+
  scale_fill_manual(values=far_colours,limits=far_levels,drop=FALSE,name='Observed FAR',
   guide=guide_legend(order=1,nrow=1,override.aes=list(fill=unname(far_colours),color=NA,alpha=1)))+
  geom_sf(data=visible,fill=NA,color='white',linewidth=1.25)+geom_sf(data=visible,aes(color=rcoi_type),fill=NA,linewidth=.65)+
  scale_color_manual(values=type_colours,limits=used,drop=FALSE,name='RCOI boundary type',guide=guide_legend(order=2,nrow=2,override.aes=list(linewidth=1.2)))+
  geom_sf(data=boundary,fill=NA,color='#555555',linetype='dashed',linewidth=.3)+geom_sf(data=windows,fill=NA,color='#111111',linetype='dashed',linewidth=.6)+
  geom_sf(data=centre,shape=23,fill='#FFD700',color='#111111',size=2.7)+
  coord_sf(crs=3035,datum=NA,xlim=ee[c('xmin','xmax')],ylim=ee[c('ymin','ymax')],expand=FALSE)+
  labs(title='(a) FAR and subcentres',subtitle=paste0('40 × 40 km | 500 m FAR | ',secondary_count,' secondary patches with at least 4 cells'))+theme_map
 for(i in 1:3){bb<-st_bbox(windows[i,]);overview<-overview+annotate('label',x=bb['xmin'],y=bb['ymax'],label=as.character(i),fontface='bold',fill='white',size=3.5,hjust=0,vjust=1)}
 overview<-adorn(overview,ee,5000)
 details<-list();points_all<-list();counts<-list()
 st_write(patches,gpkg,layer='classified_threshold_patches',delete_layer=TRUE,quiet=TRUE)
 st_write(windows,gpkg,layer='detail_windows',delete_layer=TRUE,quiet=TRUE)
 st_write(far,gpkg,layer='observed_FAR_overview_cells',delete_layer=TRUE,quiet=TRUE)
 for(i in 1:3){
  w<-windows[i,];bb<-st_bbox(w);target<-patches[patches$patch_id==w$patch_id,]
  context<-get_context(city,w)
  points<-st_read(pm$poi_path,query='SELECT osm_uid,name,commercial,employment,industrial_poi,representative_point,geom FROM all_classified_pois WHERE commercial=1 OR employment=1 OR industrial_poi=1',
   wkt_filter=st_as_text(st_geometry(w)),quiet=TRUE)
  points<-points[lengths(st_intersects(points,w))>0,];stopifnot(!anyDuplicated(points$osm_uid))
  long<-bind_rows(lapply(1:3,function(j){flag<-c('commercial','industrial_poi','employment')[j];z<-points[as.logical(points[[flag]]),];z$poi_type<-rep(poi_names[j],nrow(z));z}))
  long$detail_area<-rep(i,nrow(long));long$study_city<-rep(city,nrow(long))
  surrounding<-patches[lengths(st_intersects(patches,w))>0,]
  p<-ggplot()+geom_sf(data=context$buildings,fill='#DFE3E6',color=NA)+geom_sf(data=context$roads,color='#B8BFC4',linewidth=.15)+
   geom_sf(data=surrounding,fill=NA,color=surrounding$rcoi_colour,linewidth=.45,linetype='dashed')+
   geom_sf(data=target,fill=target$rcoi_colour,alpha=.05,color=NA)+geom_sf(data=target,fill=NA,color=target$rcoi_colour,linewidth=.95)+
   geom_sf(data=long,aes(color=poi_type,shape=poi_type,size=poi_type),alpha=.85,stroke=.35,show.legend=TRUE)+
   scale_color_manual(values=poi_colours,limits=poi_names,drop=FALSE,name='POI category')+
   scale_shape_manual(values=poi_shapes,limits=poi_names,drop=FALSE,name='POI category',guide='none')+
   scale_size_manual(values=poi_sizes,limits=poi_names,drop=FALSE,guide='none')+
   guides(color=guide_legend(nrow=1,override.aes=list(size=2,alpha=1,shape=unname(poi_shapes),colour=unname(poi_colours),stroke=.35)))+
   coord_sf(crs=3035,datum=NA,xlim=bb[c('xmin','xmax')],ylim=bb[c('ymin','ymax')],expand=FALSE)+
   labs(title=if(i==1)'(b) Commercial, industrial and office POIs'else NULL,
    subtitle=paste0(i,' · ',target$patch_label,' · ',target$rcoi_type,' | 2 × 2 km | C: ',sum(points$commercial),'  I: ',sum(points$industrial_poi),'  O: ',sum(points$employment)))+theme_map
  details[[i]]<-adorn(p,bb,500,2.3);points_all[[i]]<-long
  counts[[i]]<-data.frame(city=city,detail_area=i,patch_id=target$patch_id,patch_label=target$patch_label,cells=target$cells,rcoi_type=target$rcoi_type,
   boundary_colour=target$rcoi_colour,commercial=sum(points$commercial),industrial=sum(points$industrial_poi),office_employment=sum(points$employment),
   unique_CIO_objects=nrow(points),plotted_category_memberships=nrow(long),buildings=nrow(context$buildings),road_features=nrow(context$roads),
   detail_windows_overlap=chosen$overlap,xmin=unname(bb['xmin']),ymin=unname(bb['ymin']),xmax=unname(bb['xmax']),ymax=unname(bb['ymax']))
  st_write(context$buildings,gpkg,layer=paste0('buildings_area_',i),delete_layer=TRUE,quiet=TRUE)
  st_write(context$roads,gpkg,layer=paste0('roads_area_',i),delete_layer=TRUE,quiet=TRUE)
 }
 if(minimal_text){
  overview<-overview+labs(title=NULL,subtitle=NULL)
  details<-lapply(details,function(p)p+labs(title=NULL,subtitle=NULL))
 }
 right<-wrap_plots(details,ncol=1,guides='collect')&theme(legend.position='bottom')
 figure<-(overview|right)+plot_layout(widths=c(1.65,1))+plot_annotation(title=city,caption=if(minimal_text)NULL else paste0(
  if(above_mean)'Subcentres: FAR >= 30% of central FAR; at least 4 cells; above the city mean in at least one R/C/O/I category. Original city coverage retained. Main view: 40 × 40 km.\n'else
   'Subcentres: observed FAR >= 30% of central FAR; patches with at least 4 cells; the central core is excluded from the >5 secondary-patch city-selection count. Main view: 40 × 40 km.\n',
  'RCOI boundary types: patch mean cell POI counts standardized against city-wide category means/SDs, including zero-count cells. R: address proxy; C: commercial; O: office/employment; I: industrial.\n',
  'Right: three 2 × 2 km examples ranked by commercial POIs; symbols indicate C (circle), I (triangle), O (open square). Multi-category objects can carry several symbols. All visible patch outlines meet the four-cell rule.\n',
  'FAR: shared red-orange classes; gray means unavailable FAR. Gold diamond: CBD. Two-tone compass arrows indicate true north. Sources: saved FAR, OSM/Geofabrik POIs and roads, EU/OSM buildings, INSEE.'),
  theme=theme(plot.title=element_text(size=22,face='bold',hjust=.5),plot.caption=element_text(size=8,hjust=0),plot.margin=margin(10,10,10,10)))
 output<-file.path(out,'maps',paste0(stem,'_FAR_RCOI_CIO_min4_redorange'))
 ggsave(paste0(output,'.png'),figure,width=50,height=38,units='cm',dpi=300,device=ragg::agg_png,bg='white')
 ggsave(paste0(output,'_preview.png'),figure,width=50,height=38,units='cm',dpi=100,device=ragg::agg_png,bg='white')
 ggsave(paste0(output,'.pdf'),figure,width=50,height=38,units='cm',device=pdf,useDingbats=FALSE,encoding='WinAnsi.enc',bg='white')
 st_write(bind_rows(points_all),gpkg,layer='commercial_industrial_office_points',delete_layer=TRUE,quiet=TRUE)
 write_csv(st_drop_geometry(patches),file.path(out,'maps',paste0(stem,'_retained_RCOI_patches.csv')))
 list(plot=figure,patches=patches,excluded=excluded,reference=classified$reference,details=bind_rows(counts),summary=data.frame(city=city,retained_secondary_patches=secondary_count,
  above_mean_filter=above_mean,minimal_text=minimal_text,excluded_patches=nrow(excluded),still_more_than_5_subcentres=secondary_count>5,
  central_components=sum(patches$is_central_core),visible_secondary_patches=sum(!visible$is_central_core),unavailable_poi_patches=sum(patches$analysed_cells==0),
  partially_covered_poi_patches=sum(patches$poi_cell_coverage_fraction>0 & patches$poi_cell_coverage_fraction<1),overview_width_km=40,detail_width_km=2,
  north_arrow='Two-tone true-north compass',detail_windows_overlap=chosen$overlap,output_pdf=paste0(output,'.pdf')))
}
pdf(file.path(out,'All_eligible_cities_FAR_RCOI_CIO_redorange.pdf'),width=50/2.54,height=38/2.54,useDingbats=FALSE,encoding='WinAnsi.enc')
results<-list();status<-list()
for(city in cities){
 z<-tryCatch(make_city(city),error=function(e){message('FAILED ',city,': ',conditionMessage(e));e})
 if(inherits(z,'error'))status[[city]]<-data.frame(city=city,status='failed',detail=conditionMessage(z))else{
  print(z$plot);z$plot<-NULL;results[[city]]<-z;status[[city]]<-data.frame(city=city,status='complete',detail=z$summary$output_pdf)
 }
 write_csv(bind_rows(status),file.path(out,'run_status.csv'))
 if(length(results)){
  write_csv(bind_rows(lapply(results,`[[`,'summary')),file.path(out,'map_manifest.csv'))
  write_csv(bind_rows(lapply(results,`[[`,'details')),file.path(out,'all_detail_POI_counts.csv'))
  write_csv(bind_rows(lapply(results,`[[`,'reference')),file.path(out,'city_RCOI_references.csv'))
 }
 rm(z);gc(FALSE)
}
dev.off()
if(any(bind_rows(status)$status!='complete'))stop('Some maps failed; see run_status.csv')
st_write(bind_rows(lapply(results,`[[`,'patches')),file.path(out,'All_eligible_cities_RCOI_min4.gpkg'),layer='classified_threshold_patches',delete_layer=TRUE,quiet=TRUE)
if(above_mean)write_csv(st_drop_geometry(bind_rows(lapply(results,`[[`,'excluded'))),file.path(out,'excluded_patches.csv'))
message('COMPLETE: ',length(cities),' city figures: ',normalizePath(out))
