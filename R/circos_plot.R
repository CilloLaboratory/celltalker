#' Creates a circos plot from the list of ligands and receptors
#'
#' @param ligand_receptor_frame Resulting tibble (usually filtered in some way)
#' from the celltalk function.
#'
#' @param cell_group_colors Colors used for the groups of cells in the outer
#' track of the circos plot.
#'
#' @param ligand_color Color to use for ligands. Defaults to "blue".
#'
#' @param receptor_color Color to use for the receptors. Defaults to "red".
#'
#' @param cex_outer Size of the text for the cell groups in the outer layer of
#' the circos plot. Default is 0.5.
#'
#' @param cex_inner Size of the text for the ligand and receptors in the
#' inner layer of the circos plot. Default is 0.4.
#'
#' @param facing_outer Direction of the text in the outer layer. Options are "inside", "outside", "clockwise", "reverse.clockwise", "downward", "bending.inside" and "bending.outside". Default is "downward".
#' 
#' @param facing_inner Direction of the text in the outer layer. Options are "inside", "outside", "clockwise", "reverse.clockwise", "downward", "bending.inside" and "bending.outside". Default is "downward".
#' 
#' @param lwd Control line thickness in the circos plots. Default is 3.
#' 
#' @param arr.length Control the length of the arrow head in circos plots. Default is 0.2
#' 
#' @param arr.width Control the width of the arrow head in circos plots. Default is half of 10% of the arrow width.
#' 
#' @return Generates a circos plot connecting ligands and receptors across cell types for a given sample group
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr group_split
#' @importFrom dplyr summarize
#' @importFrom dplyr arrange
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr group_by
#' @importFrom dplyr distinct
#' @importFrom tibble enframe
#' @importFrom circlize CELL_META
#' @importFrom circlize circos.link
#' @importFrom circlize circos.track
#' @importFrom circlize circos.clear
#' @importFrom circlize circos.par
#' @importFrom circlize circos.initialize
#' @importFrom circlize circos.rect
#' @importFrom circlize circos.text
#' @importFrom circlize get.cell.meta.data
#'
#' @export

circos_plot <- function(ligand_receptor_frame,
  cell_group_colors,
  ligand_color="blue",
  receptor_color="red",
  cex_outer=0.5,
  cex_inner=0.4,
  facing_outer="downward",
  facing_inner="downward",
  lwd=3,
  arr.length=0.2, 
  arr.width=(3*0.1)/2) {

## Pieces to draw outer rec
ligand_cell_types <- ligand_receptor_frame %>%
  split(.$cell_type1) %>%
  sapply(.,nrow) %>%
  enframe(.,name="cell_type",value="lig_length")

receptor_cell_types <- ligand_receptor_frame %>%
  split(.$cell_type2) %>%
  sapply(.,nrow) %>%
  enframe(.,name="cell_type",value="rec_length")

outer_frame_length <- full_join(ligand_cell_types, receptor_cell_types,
                                by="cell_type") %>%
  mutate(sum = rowSums(across(where(is.numeric)),na.rm=T)) %>%
  split(.$cell_type) %>%
  lapply(.,function(x)
  {
    if (x$sum>1) {
    data.frame(cell_type=rep(x$cell_type, x$sum),
                x_vec=1:x$sum)
    } else {
        data.frame(cell_type=rep(x$cell_type, x$sum+1),
                x_vec=1:(x$sum+1))
    }
  }) %>%
  do.call(rbind, .)

## Pieces to draw rectangles
ligands <- ligand_receptor_frame %>%
  mutate(lig_rec=sapply(strsplit(interaction,split="_"),function(x) x[[1]])) %>%
  mutate(lig_rec_class="ligand") %>%
  select(cell_type1,lig_rec,lig_rec_class) %>%
  distinct()
  
receptors <- ligand_receptor_frame %>%
  mutate(lig_rec=sapply(strsplit(interaction,split="_"),function(x) x[[2]])) %>%
  mutate(lig_rec_class="receptor") %>%
  select(cell_type2,lig_rec,lig_rec_class) %>%
  distinct()

colnames(ligands)[1] <- colnames(receptors)[1] <- c("cell_type")

ligands_receptors <- rbind(ligands,receptors) %>% 
  group_by(cell_type) %>%
  mutate(x_vec=row_number()) %>%
  ungroup()

## Add a row for any singlets
cell_type_numbers <- ligands_receptors %>%
  split(.$cell_type) %>%
  sapply(.,function(x) nrow(x))

if (any(cell_type_numbers==1)) {
  singlet_index <- which(cell_type_numbers==1)
  for (i in 1:length(singlet_index)) {
    new_row <- data.frame(cell_type=names(singlet_index)[1],
      lig_rec=NA,
      lig_rec_class=NA,
      x_vec=2
    )
    ligands_receptors <- rbind(ligands_receptors,new_row)
  }
}

# inner_track_length <- full_join(ligands_receptors,outer_frame_length,by=c("cell_type",# "x_vec")) %>%
#  mutate(lig_rec=ifelse(is.na(lig_rec),"XXX",lig_rec)) %>%
#  mutate(lig_rec_class=ifelse(is.na(lig_rec_class),"XXX",lig_rec_class)) %>%
#  arrange(cell_type)

inner_track_length <- ligands_receptors %>%
  arrange(cell_type)

## Initalize circos plot 
inner_track_length$cell_type <- as.factor(inner_track_length$cell_type)

circos.clear()
circos.par(gap.degree = 10, track.margin=c(0,0.2))
circos.initialize(sectors=inner_track_length$cell_type, x = inner_track_length$x_vec)

## Build outer layer 
circos.track(ylim = c(0, 2.5), track.height = 0.2, bg.border=NA, panel.fun = function(x, y) {
  circos.rect(xleft=CELL_META$cell.xlim[1], xright=CELL_META$cell.xlim[2], ybottom=CELL_META$cell.ylim[1],ytop=CELL_META$cell.ylim[2]-1.5,
              col = cell_group_colors[CELL_META$sector.numeric.index])
  circos.text(x=CELL_META$xcenter,y=2.5,CELL_META$sector.index,facing = facing_outer,cex=cex_outer,niceFacing = TRUE)
})

## Build inner layer 
circos.track(ylim = c(0, 1), track.height = 0.1, bg.border=NA)
inner_track_length_split <- inner_track_length %>%
  filter(!lig_rec==is.na(NA)) %>%
  split(.$cell_type)

suppressMessages({
for (i in 1:length(inner_track_length_split)) {

  cell_type_lig_rec <- inner_track_length_split[[i]]
  sector_xlim <- get.cell.meta.data("xlim",sector.index=cell_type_lig_rec$cell_type[1])
  sector_factor <- (sector_xlim[2] - sector_xlim[1]) / nrow(cell_type_lig_rec)

  for (a in 1:nrow(cell_type_lig_rec)) {
    if (a==1) {
      circos.rect(xleft=CELL_META$xlim[1],xright=CELL_META$xlim[1]+a*sector_factor,ybottom=0,ytop=1,sector.index=cell_type_lig_rec$cell_type[a],
      col=ifelse(cell_type_lig_rec$lig_rec_class[a]=="ligand",ligand_color,receptor_color))
      circos.text(x=CELL_META$xlim[1]+a*sector_factor-(sector_factor*0.5),y=2.5,label=cell_type_lig_rec$lig_rec[a],sector.index=cell_type_lig_rec$cell_type[a],facing = facing_inner,cex=cex_inner,niceFacing = TRUE)
    } else {
      circos.rect(xleft=CELL_META$xlim[1]+(a-1)*sector_factor,xright=CELL_META$xlim[1]+a*sector_factor,ybottom=0,ytop=1,sector.index=cell_type_lig_rec$cell_type[a],
      col=ifelse(cell_type_lig_rec$lig_rec_class[a]=="ligand",ligand_color,receptor_color))
      circos.text(x=CELL_META$xlim[1]+(a)*sector_factor-(sector_factor*0.5),y=2.5,label=cell_type_lig_rec$lig_rec[a],sector.index=cell_type_lig_rec$cell_type[a],facing = facing_inner,cex=cex_inner,niceFacing = TRUE)
    }
  }
}
})

## Build links
link_frame <- ligand_receptor_frame %>% 
  mutate(ligand=sapply(strsplit(interaction,split="_"),function(x) x[[1]])) %>%
  mutate(receptor=sapply(strsplit(interaction,split="_"),function(x) x[[2]])) %>%
  select(cell_type1,cell_type2,ligand,receptor)

for (i in 1:nrow(link_frame)) {

test <- link_frame[i,]

lig_index <- ligands_receptors %>%
  filter(cell_type %in% test$cell_type1) %>%
  filter(lig_rec %in% test$ligand) %>%
  filter(lig_rec_class=="ligand") %>%
  pull(x_vec)

cell_type_lig <- ligands_receptors %>%
  filter(!is.na(lig_rec)) %>%
  filter(cell_type %in% test$cell_type1)

sector_xlim_lig <- get.cell.meta.data("xlim",sector.index=test$cell_type1)
sector_factor_lig <- (sector_xlim_lig[2] - sector_xlim_lig[1]) / nrow(cell_type_lig)
lig_loc <- sector_xlim_lig[1]+(lig_index)*sector_factor_lig-(sector_factor_lig*0.5)

rec_index <- ligands_receptors %>%
  filter(cell_type %in% test$cell_type2) %>%
  filter(lig_rec %in% test$receptor) %>%
  filter(lig_rec_class=="receptor") %>%
  pull(x_vec)

cell_type_rec <- ligands_receptors %>%
  filter(!is.na(lig_rec)) %>%
  filter(cell_type %in% test$cell_type2)

sector_xlim_rec <- get.cell.meta.data("xlim",sector.index=test$cell_type2)
sector_factor_rec <- (sector_xlim_rec[2] - sector_xlim_rec[1]) / nrow(cell_type_rec)
rec_loc <- sector_xlim_rec[1]+(rec_index)*sector_factor_rec-(sector_factor_rec*0.5)

suppressMessages({
  circos.link(sector.index1 = test$cell_type1, point1=lig_loc, sector.index2 = test$cell_type2, point2=rec_loc, rou1=0.3, rou2=0.3, directional=1, lwd=3, arr.length=0.2, arr.width=(3*0.1)/2)
})
}

}