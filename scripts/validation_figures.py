"""
==============================================================================
Functions to reproduce figures from Paquin et al. (Scientific Data)
Author: Christopher McCray (Ouranos)
Date: April 2025
==============================================================================
This file contains a set of functions used to produce the figures 
presented in Paquin et al. (Scientific Data). It includes tools for generating
multi-panel seasonal maps, annual cycles, boxplots, and region-specific 
comparative visualizations of climate simulation datasets.

- `postage_fig`:
    Generate subplots using rotated pole projection with customizable extent.

- `create_figure_layout`:
    Helper function to generate layouts for common figure orientations.

- `affiche_cycle_annuel`:
    Plot annual cycles and differences between two datasets for a given region.

- `cartes_differences_sim_ref_8panels`:
    Generate an 8-panel figure comparing seasonal differences for two variables.

- `boxplots_saisons`:
    Create seasonal boxplots or violin plots comparing multiple simulations.

- `get_variable_props`:
    Determine color scales and unit conversions for plotting variables.

- `get_past_fut_mrcc5_pilote`:
    Retrieve simulation datasets from a project catalog, including past and 
    future runs under different SSPs.

- `cos_weight_mean`:
    Calculate cosine-latitude-weighted averages over specified dimensions.

- `colors_ipcc` and `pcolor_map`:
    IPCC-compliant color schemes and utilities for colormap generation.
"""
import matplotlib as mpl
import matplotlib.colors as mplc
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import numpy as np
from clisops.core import subset
import xscen as xs
import xclim as xc
import xarray as xr
import seaborn as sns
import pandas as pd
import os


def postage_fig(nrows, ncols, extent=None, figsize=(8,12), hspace=-0.5, ds_lakes=None, latlon=False, mask_lakes=True, lake_color = 'white'):
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
   
    proj = ccrs.RotatedPole(pole_longitude=83, pole_latitude=42.5)  
    if extent == None:
        #extent = [-128,-62,10.5,83]   
        extent = [-34.04499817, 37.89500427, -33.62499619, 35.34500504]  
        extent = [-34.2, 38, -33.7, 35.5] 
    zoom_proj = proj   
    #zoom_proj = ccrs.PlateCarree(central_longitude=0)
    #proj.threshold = 500
    #zoom_proj.threshold=500
    
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,  subplot_kw=dict(projection=proj))
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face',
                                            facecolor='white')   
    
    for ax in axs:
        if len(axs.shape)==1:        
            naxes = 1
            ax = [ax]
        else:
            naxes = len(ax)
        for nax in np.arange(0,naxes):
            for spine in ax[nax].spines.values():
                spine.set_zorder(999)
            ax[nax].set_extent(extent,zoom_proj)  
            #ax[1].set_extent(extent,ccrs.PlateCarree())  
            ax[nax].coastlines('50m', zorder=10, linewidths=0.5)
            
            states_provinces = cfeature.NaturalEarthFeature(
                category='cultural',  name='admin_1_states_provinces_lines',
                scale='50m', facecolor='none')

            ax[nax].add_feature(land_50m,facecolor='0.9')            
            ax[nax].add_feature(cfeature.BORDERS, zorder=10, linewidths=0.5)
            ax[nax].add_feature(states_provinces, edgecolor='0.5', zorder=9, linewidths=0.5)
            ax[nax].add_feature(cfeature.LAKES,edgecolor='k', facecolor='none', zorder=9, linewidths=0.5)
            
            if ds_lakes is not None:
                lf_thresh = 95 
                lakefrac = ds_lakes.where(ds_lakes.lakeFrac>=lf_thresh).lakeFrac
                if mask_lakes == True:
                    ax[nax].contourf(ds_lakes.lon, ds_lakes.lat, lakefrac, transform=ccrs.PlateCarree(), zorder=10, colors=lake_color)
                ax[nax].contour(ds_lakes.lon, ds_lakes.lat, ds_lakes.lakeFrac.values, transform=ccrs.PlateCarree(),                    zorder=11, colors='k', linewidths=0.5, levels=[lf_thresh])
            
            #Next 3 lines, lat/lon lines
            if latlon==True:
                gl = ax[nax].gridlines(crs=zoom_proj, linewidth=1, color='black', alpha=0.5, linestyle='--',zorder=99)
                gl.xlocator = mticker.FixedLocator(np.arange(-180,181,20))
                gl.ylocator = mticker.FixedLocator(np.arange(0,91,20))

    fig.subplots_adjust(wspace=0.05, hspace=hspace)
    return fig, axs

def create_figure_layout(orientation='2x2'):
    '''
    Creates a figure and axes based on the given orientation.
    '''
    if orientation == '2x2':
        fig, axs = postage_fig(2, 2, figsize=(12, 16), hspace=-0.35)
        cbar_axes = [0.92, 0.31, 0.02, 0.39]
        y_title = 0.85
        fontsize=22
    elif orientation == '1x4':
        fig, axs = postage_fig(1, 4, figsize=(19, 7), hspace=-0.35)
        cbar_axes = [0.92, 0.25, 0.012, 0.49]
        y_title = 0.90
        fontsize=20
    elif orientation == '2x4':
        fig, axs = postage_fig(2, 4, figsize=(19, 14), hspace=-0.45)
        cbar_axes = [0.92, 0.31, 0.02, 0.39]
        y_title = 0.85
        fontsize=22
    else:
        raise ValueError("Invalid orientation. Use '2x2','1x4' or '2x4'.")
    return fig, axs, cbar_axes, y_title, fontsize


def affiche_cycle_annuel(ds_ref, 
                         ds_source, 
                         region, 
                         nom_region, 
                         nom_var, 
                        fr_en = 'fr',
                        label_attr = 'code'):
    lw = 3
    fontsize = 14
    def delta_colors(nom_var):
        if nom_var.startswith('tasmax'):
            delta_col = 'tab:red'
        elif nom_var.startswith('tasmin'):
            delta_col = 'tab:blue'
        elif nom_var.startswith('tas'):
            delta_col = '0.1'
        return delta_col
    
    ds_ref = subset.subset_shape(ds_ref, shape=region)
    ds_source = subset.subset_shape(ds_source, shape=region)

    ds_reg_ref = xs.spatial_mean(ds_ref, method='cos-lat')
    ds_reg_source = xs.spatial_mean(ds_source, method='cos-lat')
    
    if 'experiment_id' in ds_reg_ref.attrs:
        if label_attr == 'code':
            nom_ref = ds_reg_ref.attrs['experiment_id']    
        else:
            if 'driving_model_id' in ds_reg_ref.attrs:
                pilote = ds_reg_ref.attrs['driving_model_id'].strip()
            else:
                pilote = ds_reg_ref.attrs['cat:driving_model']
            nom_ref = f"CRCM5[{pilote}]"
    else:
        nom_ref = ds_reg_ref.attrs['cat:source']
        
    if 'experiment_id' in ds_reg_source.attrs:
        if label_attr == 'code':
            nom_source = ds_reg_source.attrs['experiment_id']      
        else:
            if 'driving_model_id' in ds_reg_source.attrs:
                pilote = ds_reg_source.attrs['driving_model_id'].strip()
            else:
                pilote = ds_reg_source.attrs['cat:driving_model']
            nom_source = f"CRCM5[{pilote}]"         
    else:
        nom_source = ds_reg_source.attrs['cat:source']
        
    periode_ref = ds_reg_ref[f"{nom_var}_moy_mens"].period
    periode_source = ds_reg_source[f"{nom_var}_moy_mens"].period

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=[5,8])
    ax_values = axs.flatten()[0]
    ax_delta = axs.flatten()[1]
    #Valeurs générales            
    mois = ds_reg_ref.month.values
    #nom_reg = ds_reg_ref.Name.values    
    colors = []
    # Faire graphique avec les valeurs de ref et source
    for ref_source, ds_reg in zip(['ref','source'], [ds_reg_ref, ds_reg_source]):
        if 'cat:experiment' in ds_reg.attrs:
            exp = ds_reg.attrs['cat:experiment']
        else:
            exp = None

        c = colors_ipcc(exp)
        # S'assurer que les 2 couleurs ne sont pas les mêmes
        colors.append(c)
        if len(colors) > 1:
            if colors[0] == colors[1]:
                c = '0.3'
        
        if 'experiment_id' in ds_reg.attrs:
            #source_label = ds_reg.attrs['experiment_id']
            source_label = ds_reg.attrs['experiment_id']
        else:
            source_label = ds_reg.attrs['cat:source']
        period = ds_reg[f"{nom_var}_moy_mens"].period
        leg_lab = f"{source_label} {period}"
        
        if nom_var.startswith('tas'):
            tas_moy = ds_reg.tas_moy_mens.values-273.15
            ax_values.plot(mois,  tas_moy, label=leg_lab, c=c, lw=lw)
            if 'tasmax_moy_mens' in ds_reg.data_vars:
                tasmax_moy = ds_reg.tasmax_moy_mens.values-273.15
                tasmin_moy = ds_reg.tasmin_moy_mens.values-273.15            
                leg_lab_fill = f"tasmax-tasmin"
                ax_values.plot(mois,  tasmax_moy, label=None, c=c, lw=1, ls='--')
                ax_values.plot(mois,  tasmin_moy, label=None, c=c, lw=1, ls='--')
                ax_values.fill_between(mois, tasmin_moy, tasmax_moy, alpha=0.1, color=c, 
            label=leg_lab_fill)                
                tas_vars = ['tas', 'tasmin','tasmax']
            else:
                tas_vars = ['tas']
            #plt.axhline(y=0, c= 'k')
            ylabel = '˚C'
        elif nom_var.startswith('prsn'):
            prsn_moy = ds_reg.prsn_moy_mens.values * 86400
            ax_values.plot(mois,  prsn_moy, label=leg_lab, c=c, lw=lw)
            #plt.ylim(0,14)
            ylabel = r'mm j$^{-1}$'
        elif nom_var.startswith('pr'):
            pr_moy = ds_reg.pr_moy_mens.values * 86400
            ax_values.plot(mois,  pr_moy, label=leg_lab, c=c, lw=lw)
            #plt.ylim(0,14)
            if fr_en == 'fr':
                ylabel = r'mm j$^{-1}$'
            elif fr_en == 'en':
                ylabel = r'mm day$^{-1}$'
    if fr_en == 'fr':
        valeurs_title = f"{nom_region}\nCycle annuel - {nom_var}"     
        diff = 'Différence'
    elif fr_en == 'en':
        valeurs_title = f"{nom_region}\nSeasonal cycle - {nom_var}"
        diff = 'Difference'
    if periode_source == periode_ref:
        deltas_title = f"{diff}\n{nom_source} - {nom_ref} ({periode_ref})"
    else:
        deltas_title = f"{diff}\n{nom_source} ({periode_source}) - {nom_ref} ({periode_ref})"
    ax_values.set_title(valeurs_title, fontsize=fontsize)
    ax_values.legend()


    # Faire graphique avec les deltas

    if nom_var.startswith('prsn'):
        leg_lab = fr'$\Delta${nom_var}'
        delta =(( ds_reg_source[f'{nom_var}_moy_mens'] - ds_reg_ref[f'{nom_var}_moy_mens']) * 86400 ).values
        #ylim = [-3,3]
        ax_delta.plot(mois, delta, c='0.1', lw=lw, alpha=0.9, label=leg_lab)
    elif nom_var.startswith('pr'):
        leg_lab = fr'$\Delta${nom_var}'
        delta =(( ds_reg_source[f'{nom_var}_moy_mens'] - ds_reg_ref[f'{nom_var}_moy_mens']) * 86400 ).values
        #ylim = [-3,3]
        ax_delta.plot(mois, delta, c='0.1', lw=lw, alpha=0.9, label=leg_lab)
    else:
        #tas
        #ylim = [-10,10]
        for tas_var in tas_vars:
            leg_lab = fr'$\Delta${tas_var}'
            delta =( ds_reg_source[f'{tas_var}_moy_mens'] - ds_reg_ref[f'{tas_var}_moy_mens']).values
            ax_delta.plot(mois, delta, c=delta_colors(tas_var), lw=lw, alpha=0.9, label=leg_lab)            

    ax_delta.axhline(y=0, c= 'k')
   
    ax_delta.set_title(deltas_title, fontsize=fontsize)
    ax_delta.legend()



    for ax in axs.flatten():            
        ax.grid()              
        ax.set_xlim([1,12])
        ax.set_xticks(np.arange(1,13))
        #ax.set_ylabel(ylabel, fontsize=fontsize-2, rotation=0, labelpad=10)
        ax.set_ylabel(ylabel, fontsize=fontsize-2, )
        ax.tick_params(axis='both', which='major', labelsize=fontsize-3)
        #ax.set_xticks(ax.get_xticks())
        #ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        #ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
        #ax.show()
    #fig.tight_layout()
    plt.subplots_adjust(hspace=0.3)
    
    #left, bottom, width, height = [0.77, 0.9, 0.125, 0.125]
    left, bottom, width, height = [0.77, 0.87, 0.125, 0.125]
    proj = ccrs.RotatedPole(pole_longitude=83, pole_latitude=42.5)  
    ax2 = fig.add_axes([left, bottom, width, height], projection=proj)
    proj = ccrs.RotatedPole(pole_longitude=83, pole_latitude=42.5)  
    ax2.set_extent([-128,-62,10.5,83],crs=ccrs.PlateCarree())  
    #states_provinces = cfeature.NaturalEarthFeature(
    #    category='cultural',  name='admin_1_states_provinces_lines',
    #    scale='50m', facecolor='none')
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face',
                                            facecolor='white')
    ax2.add_feature(land_50m);
    ax2.add_feature(cfeature.BORDERS, zorder=10); 
    #ax2.add_feature(states_provinces, edgecolor='0.5', zorder=9)
    ax2.coastlines('110m', zorder=10)
    ax2.add_feature(cfeature.LAKES,edgecolor='k', facecolor='none', zorder=9, linewidths=0.5)
    
    #region.plot(ax=ax2, transform=ccrs.PlateCarree(), zorder=5)
    plot_data =  ds_ref[f"{nom_var}_moy_mens"].isel(month=0).notnull()
    cmap = mplc.ListedColormap(['None','purple'])
    cb = ax2.pcolormesh(plot_data.lon, plot_data.lat, plot_data*1,vmin=0.5,vmax=1.1,
                 zorder=5, transform=ccrs.PlateCarree(), cmap=cmap)
    return fig
  
        
def cartes_differences_sim_ref_8panels(comp, 
                                      ds_ref, 
                                      ds_source, 
                                      nom_var1, 
                                      nom_var2,  # New parameter for second variable
                                      abs_rel_pr = 'rel', 
                                      ref_period=None, 
                                      source_period=None, 
                                      levels=None, 
                                      max_min=True,
                                      title=True):
    xr.set_options(keep_attrs=True)
    
    fig, axs, cbar_axes, y_title, fontsize = create_figure_layout('2x4')  # Set to 2x4 layout for 8 panels

    # Extract the first part of the variable names for labeling on the left
    if 'ERA5' in ds_source.attrs['driving_model_id']:
        pilote = 'ERA5'
    else:
        pilote = ds_source.attrs['driving_model_id']
    text_diff =  f" (CRCM5[{pilote}]–{ds_ref.attrs['cat:source']})"
    nom_var1_label = nom_var1.split('_')[0]+text_diff
    nom_var2_label = nom_var2.split('_')[0]+text_diff

    # Define labels for the top and bottom rows
    labels_top_row = ['a', 'b', 'c', 'd']
    labels_bottom_row = ['e', 'f', 'g', 'h']
    
    # First variable (top row)
    da_var_ref1 = ds_ref[nom_var1]
    da_var_source1 = ds_source[nom_var1]

    # Second variable (bottom row)
    da_var_ref2 = ds_ref[nom_var2]
    da_var_source2 = ds_source[nom_var2]
    abs_rel_var1 = 'abs'
    abs_rel_var2 = 'abs'
    da_var1 = da_var_source1 - da_var_ref1
    da_var2 = da_var_source2 - da_var_ref2
    
    if 'pr' in nom_var1:
        abs_rel_var1 = abs_rel_pr
        if abs_rel_pr  == 'rel':
            da_var1 = (da_var_source1 - da_var_ref1) / da_var_ref1 * 100
        else:
            da_var1 = da_var_source1 - da_var_ref1
    if 'pr' in nom_var2:
        abs_rel_var2 = abs_rel_pr
        if abs_rel_pr  == 'rel':
            da_var2 = (da_var_source2 - da_var_ref2) / da_var_ref2 * 100
        else:
            da_var2 = da_var_source2 - da_var_ref2
    
    # Loop over seasons for both variables (top and bottom row)
    for n, season in enumerate(['DJF', 'MAM', 'JJA', 'SON']):
        # Top row: nom_var1 (tas_moy_sais)
        plot_ax1 = axs.flatten()[n] 
        for spine in plot_ax1.spines.values():
            spine.set_zorder(999)
                
        if season in da_var1.season:
            da_diff_season1 = da_var1.sel(season=season)
            da_diff_season1, cmap1, bounds1, norm1, extend1, units1 = get_variable_props(nom_var=nom_var1, 
                                                                                   da_var=da_diff_season1, 
                                                                                   vals_diff='diff', 
                                                                                   abs_rel=abs_rel_var1,
                                                                                        fr_en='en')
            lon = da_diff_season1.lon
            lat = da_diff_season1.lat
            cf1 = plot_ax1.pcolormesh(lon, lat, da_diff_season1, cmap=cmap1, norm=norm1, transform=ccrs.PlateCarree(), zorder=4)
            
            # Only add season titles for the top row
            if n == 0:
                plot_ax1.text(-0.1, 0.5, f"{nom_var1_label}", transform=plot_ax1.transAxes, 
                              fontsize=fontsize-2, va='center', ha='center', rotation='vertical',zorder=999)
            plot_ax1.set_title(f"{season}", fontsize=fontsize)
            
            # Add label to the top row (A-D)
            letter_label = "(" + labels_top_row[n] + ")"
            plot_ax1.text(0.04, 0.04, letter_label, transform=plot_ax1.transAxes,
                         fontsize=fontsize, ha='left', va='bottom', zorder=999,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.5, pad=2))
        
        # Bottom row: nom_var2 (pr_moy_sais)
        plot_ax2 = axs.flatten()[n + 4]  # Bottom row starts at index 4
        for spine in plot_ax2.spines.values():
            spine.set_zorder(999)
        if season in da_var2.season:
            da_diff_season2 = da_var2.sel(season=season)
            da_diff_season2, cmap2, bounds2, norm2, extend2, units2 = get_variable_props(nom_var=nom_var2, 
                                                                                   da_var=da_diff_season2, 
                                                                                   vals_diff='diff', 
                                                                                   abs_rel=abs_rel_var2,
                                                                                        fr_en='en')
            lon = da_diff_season2.lon
            lat = da_diff_season2.lat
            cf2 = plot_ax2.pcolormesh(lon, lat, da_diff_season2, cmap=cmap2, norm=norm2, transform=ccrs.PlateCarree(), zorder=4)
            
            # Add variable name on the left side of the bottom row
            if n == 0:
                plot_ax2.text(-0.1, 0.5, f"{nom_var2_label}", transform=plot_ax2.transAxes, 
                              fontsize=fontsize-2, va='center', ha='center', rotation='vertical', 
                              zorder=999)
            
            # Add label to the bottom row (E-H)
            #plot_ax2.text(0.02, 0.98, labels_bottom_row[n], transform=plot_ax2.transAxes, 
            #              fontsize=fontsize+4, va='top', ha='left')
            letter_label = "(" + labels_bottom_row[n] + ")"
            plot_ax2.text(0.04, 0.04, letter_label, transform=plot_ax2.transAxes,
                         fontsize=fontsize, ha='left', va='bottom', zorder=999,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.5, pad=2))


    # Add a vertical colorbar for the top row (tas_moy_sais)
    cbar_ax1 = fig.add_axes([0.91, 0.51, 0.012, 0.25])  # Aligned exactly with the height of the top row
    cb1 = fig.colorbar(cf1, cax=cbar_ax1, orientation='vertical', extend=extend1, ticks=bounds1)
    cb1.ax.tick_params(labelsize=fontsize-8)
    cb1.set_label(units1, fontsize=fontsize-6)
    
    # Add a vertical colorbar for the bottom row (pr_moy_sais)
    cbar_ax2 = fig.add_axes([0.91, 0.235, 0.012, 0.25])  # Aligned exactly with the height of the bottom row
    cb2 = fig.colorbar(cf2, cax=cbar_ax2, orientation='vertical', extend=extend2, ticks=bounds2)
    cb2.ax.tick_params(labelsize=fontsize-8)
    cb2.set_label(units2, fontsize=fontsize-6)
    
    if title:
        plt.suptitle(f"{nom_var1} (top) and {nom_var2} (bottom)", y=y_title, fontsize=fontsize+2)
    
    plt.show()
    return fig, axs

def boxplots_saisons(nom_var, 
                     ds_list, 
                     region, 
                     region_name,                      
                     box_ou_violin = 'box',
                    fr_en = 'fr',
                    label_attr = 'code'):
    sns.set_style("darkgrid")
    #colors = ['#984ea3','#e41a1c','#377eb8','#4daf4a',]      
    colors = []
        
    #fig, axs = plt.subplots(2,2, figsize=[8.5,7])
    fig, axs = plt.subplots(2,2, figsize=[11.1,8])

    n_ax = 0

    for season, ax in zip(['DJF','MAM','JJA','SON'], axs.flatten()):    
        #plt.figure(figsize=[4,4])   

        df_boxplots = pd.DataFrame(index=np.arange(0,200000))
        #comp_list = ['sim','obs','pilote','sim_b']    
        for ds in ds_list:
            if nom_var.startswith('tas'):
                ds[nom_var] =  xc.core.units.convert_units_to(ds[nom_var], "degC")
                units = ds[nom_var].attrs['units']
                if fr_en == 'fr':
                    nom_long_var = ds[nom_var].attrs['long_name']
                elif fr_en == 'en':
                    nom_long_var = nom_var.split('_')[0]+ " (seasonal mean)"
                    
            elif nom_var.startswith('pr'):
                # ds[nom_var] = ds[nom_var].copy(data=ds[nom_var].values*86400)
                #ds[nom_var].attrs['units'] = r'mm jour$^{-1}$'
                if fr_en == 'fr':
                    units = r'mm jour$^{-1}$'
                    nom_long_var = ds[nom_var].attrs['long_name']
                elif fr_en == 'en':
                    units = r'mm day$^{-1}$'
                    nom_long_var = nom_var.split('_')[0] + " (seasonal mean)"

            #mask = create_mask(x_dim=ds.lon, y_dim=ds.lat, poly=polys)
            ds_region = subset.subset_shape(ds[nom_var], shape=region)            
            values = ds_region.sel(season=season).values    
            if nom_var.startswith('pr'):
                values=values*86400
            period = ds_region.attrs['period']

            if 'experiment_id' in ds.attrs:
                if label_attr == 'code':
                    ds_name = f"{ds.attrs['experiment_id']}\n{period}"
                elif label_attr == 'experiment':
                    ds_name = f"{ds.attrs['cat:experiment']}\n{period}"
                elif label_attr == 'member':
                    ds_name = f"{ds.attrs['experiment_id']}({ds.attrs['cat:member']})\n{period}"
            else:
                ds_name = f"{ds.attrs['cat:source']}\n{period}"

            df_boxplots[ds_name] = pd.DataFrame(values[~np.isnan(values)])
            
            if 'cat:experiment' in ds.attrs:
                exp = ds.attrs['cat:experiment']
            else:
                exp = None
                
            colors.append(colors_ipcc(exp))
        if box_ou_violin == 'box':
            sns.boxplot(data=df_boxplots, 
                        whis=[5,95],
                        ax=ax, 
                        palette=colors)
            ax.tick_params(axis='x', labelsize=10)
        elif box_ou_violin == 'violin':
            sns.violinplot(data=df_boxplots, ax=ax, palette=colors)
        ax.set_title(season, fontsize=14)
        if n_ax in [0,2]:        
            ax.set_ylabel(units)
        n_ax+=1
        
    plt.subplots_adjust(hspace=0.3)
    left, bottom, width, height = [0.77, 0.9, 0.125, 0.125]
    proj = ccrs.RotatedPole(pole_longitude=83, pole_latitude=42.5)  
    ax2 = fig.add_axes([left, bottom, width, height], projection=proj)
    proj = ccrs.RotatedPole(pole_longitude=83, pole_latitude=42.5)  
    ax2.set_extent([-128,-62,10.5,83],crs=ccrs.PlateCarree())  
    #states_provinces = cfeature.NaturalEarthFeature(
    #    category='cultural',  name='admin_1_states_provinces_lines',
    #    scale='50m', facecolor='none')
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face',
                                            facecolor='white')
    ax2.add_feature(land_50m);
    ax2.add_feature(cfeature.BORDERS, zorder=10); 
    #ax2.add_feature(states_provinces, edgecolor='0.5', zorder=9)
    ax2.coastlines('110m', zorder=10)
    ax2.add_feature(cfeature.LAKES,edgecolor='k', facecolor='none', zorder=9, linewidths=0.5)
    plot_data =  ds_region.sel(season=season).notnull()
    cmap = mplc.ListedColormap(['None','purple'])
    cb = ax2.pcolormesh(plot_data.lon, plot_data.lat, plot_data*1,vmin=0.5,vmax=1.1,
                  zorder=5, transform=ccrs.PlateCarree(), cmap=cmap)
    if 'CRCM5' in ds.attrs['cat:source']:
        sim_title = f"CRCM5[{ds.attrs['driving_model_id'].strip()}]"
    else:
        sim_title = f"{ds.attrs['cat:source']}[{ds.attrs['cat:driving_model']}]"
    plt.suptitle(f"{nom_long_var}\n{sim_title}\nRegion: {region_name}",y=1.02, fontsize=14)
    #plt.tight_layout()
    plt.subplots_adjust(hspace=0.3, wspace=0.1) 
    plt.show()
    return fig, axs

def get_variable_props(nom_var, 
                       da_var, 
                       vals_diff, 
                       abs_rel=None,
                       fr_en='fr'):
    # Accepte une DataArray da_var, 
    # le nom de la variable nom_var, 
    # et si la différence est absolue ou relative abs_rel
    # Sortir propriétés pour les cartes 
    # Returns: cmap, bounds, norm, extend, units
    assert vals_diff in ['vals', 'diff'], f"vals_diff doit être 'vals' ou 'diff'"
    assert fr_en in ['fr', 'en'], f"fr_en doit être 'fr' ou 'en'"
    if nom_var.startswith('pr'):
        if vals_diff == 'vals':    
            cmap = 'BuGn'
            extend = 'both'   
            #Convert from kg/m2/s to mm/day
            da_var = da_var * 86400
            bounds = np.array([0,0.25,1,2.5,5,7.5,10,15])
            if fr_en == 'fr':
                units = r'mm j$^{-1}$'
            else:
                units = r'mm d$^{-1}$'
        elif vals_diff == 'diff':
            cmap = 'BrBG'
            extend = 'both'   
            if abs_rel == 'abs':
                #Convert from kg/m2/s to mm/day
                da_var = da_var * 86400
                bounds = np.array([0.1, 0.25, 0.5, 0.75, 1, 1.5, 2])
                if fr_en == 'fr':
                    units = r'mm j$^{-1}$'
                else:
                    units = r'mm d$^{-1}$'
            elif abs_rel == 'rel':    
                #bounds = np.array([-100,-50,-25,-15,-10,-5,0,5,10,15,25,50,100,150,200])
                bounds = np.array([1,5,10,15,25,50,100,150,200])
                units = '%'
                extend = 'both'
    elif nom_var.startswith('tas'):
        if vals_diff == 'vals':
            cmap = 'RdBu_r'
            extend = 'both'
            da_var = da_var-273.15
            bounds=np.arange(-35,36,5)
            bounds = np.append(bounds*-1,bounds)#.sort()
            bounds.sort()
            units = '°C'   
        elif vals_diff == 'diff':
            cmap = 'RdBu_r'
            if abs_rel == 'abs':    
                extend = 'both'          
                bounds=np.array([0.5,1,2.5,5,10,15])  
                bounds=np.array([0.5,1,2,3,4,5,10])
                units = '°C'
            elif abs_rel == 'rel':
                bounds = np.array([1,5,10,15,25,50,100,150,200])
                units = '%'
                extend = 'both' 
    else:
        # Identifier automatiquement les limites
        if vals_diff == 'vals':
            extend = 'both'  
            max_bound = np.percentile(da_var.values, 0.9)
            min_bound = np.percentile(da_var.values, 0.1)
            max_bound = max(abs(min_bound),abs(max_bound))
            step = max_bound/6
            bounds = np.arange(step, max_bound, step+step/1000000)
            cmap = 'RdBu_r'
            units = da_var.units
        elif vals_diff == 'diff':            
            if abs_rel == 'rel':
                cmap = 'RdBu_r'
                extend = 'both'
                bounds = np.array([1,5,10,15,25,50,100,150,200])
                units = '%'
            else:
                extend = 'both'  
                max_bound = np.percentile(da_var.values, 0.9)
                min_bound = np.percentile(da_var.values, 0.1)
                max_bound = max(abs(min_bound),abs(max_bound))
                step = max_bound/6
                bounds = np.arange(step, max_bound, step+step/1000000)
                cmap = 'RdBu_r'
                units = da_var.units
    if vals_diff == 'diff':        
        bounds = np.append(bounds*-1,bounds)#.sort()
        bounds.sort()
    cmap, norm = pcolor_map(bounds, cmap, extend)
    cmap = mpl.colormaps.get_cmap(cmap).copy()
    cmap.set_over("m")
    cmap.set_under("c")    
    return(da_var, cmap, bounds, norm, extend, units)

def get_past_fut_mrcc5_pilote(pilote, 
                       periode_past=[1971,2000],
                       periode_fut=[2071,2100],
                       pcat_fi='validation-crcm5-catalog.json', 
                       mrcc5_ou_pilote = 'MRCC5',
                       region = None,
                       output_ssps = ['ssp126','ssp370'],
                       member = None,
                       regrid_pilote = True,
                       obs = None,
                       months = None,
                             ):
    mrcc5_ou_pilote = mrcc5_ou_pilote.lower()
    assert mrcc5_ou_pilote in ['mrcc5', 'pilote'], f"mrcc5_ou_pilote doit être 'MRCC5' ou 'pilote'"
    
    ds_out = []
    date_start = f'{periode_past[0]}-01-01'
    date_end = f'{periode_past[-1]}-12-31'
    pcat = xs.ProjectCatalog(pcat_fi)
    if mrcc5_ou_pilote == 'mrcc5':
        processing_level = 'moy-mens-sais'
    elif mrcc5_ou_pilote == 'pilote':
        if regrid_pilote:
            processing_level = 'moy-mens-sais-regrid'
        else:
            processing_level = 'moy-mens-sais'
    if member == None:
        if 'MPI-ESM1-2-LR' in pilote:
            member = 'r1i1p1f1'
        else:
            member = ".*"
    
    # Past: select first ID because could be several with SSPs
    if mrcc5_ou_pilote == 'mrcc5':
        #source = 'OURANOS-CRCM5'        
        past_id = pcat.search(driving_model=f".*{pilote}.*", 
                              member = member,
                                    processing_level = processing_level,
                                    date_start=date_start,
                                    date_end = date_end)    
    elif mrcc5_ou_pilote == 'pilote':
        source = pilote        
        past_id = pcat.search(source=pilote, 
                                processing_level = processing_level,
                                member=member,
                                date_start=date_start,
                                date_end = date_end)
    if past_id:
        past_id = past_id.keys()[0]
    else:
        return [None]
    ds_sim_past = pcat[past_id].to_dask().load()
    
    
    if obs is not None:
        ds_obs = pcat.search(source=obs, 
                            domain = ds_sim_past.attrs['cat:domain'],
                            date_start=date_start,
                            date_end = date_end).to_dask().load()
        print("Obs trouvée, ajoutée à la liste")
        if region is not None:
            ds_obs = subset.subset_shape(ds_obs, shape=region) 
        ds_out.append(ds_obs)
        

    
    if region is not None:
        ds_sim_past = subset.subset_shape(ds_sim_past, shape=region) 
    ds_out.append(ds_sim_past)    
    
    # Simulations futures -----------------------------------------------------     
    # *** SSPs ***
    if output_ssps is not None:
        if isinstance(periode_fut, list):
            id_fut = f"stat_mens_sais_{int(periode_fut[0])}-{int(periode_fut[1])}$"
        elif isinstance(periode_fut, (float, int)):
            gwl_baseline = [1850,1900]
            id_fut = f"warminglevel-{'{:.1f}'.format(periode_fut)}vs{int(gwl_baseline[0])}-{int(gwl_baseline[1])}$"
            
        for ssp in output_ssps:   
            if mrcc5_ou_pilote == 'mrcc5':
                if pcat.search(id = id_fut, 
                                driving_model = pilote, 
                               experiment = ssp,
                               member = member,
                               processing_level = processing_level):
                    print(f'--SIMULATION {ssp} - TROUVÉE---')
                    ds_sim_fut =pcat.search(id = id_fut, 
                                            driving_model = pilote, 
                                            member = member,
                                            experiment = ssp,
                                             processing_level = processing_level).to_dask().load()
                    if region is not None:
                        ds_sim_fut = subset.subset_shape(ds_sim_fut, shape=region) 
                    ds_out.append(ds_sim_fut)
                    
                else:
                    print(f'--SIMULATION {ssp} - DONNÉES MANQUANTES---')
                    ds_out.append(None)
            elif mrcc5_ou_pilote == 'pilote':
                if pcat.search(id = id_fut, 
                                source = pilote, 
                               member = member,
                               experiment = ssp,
                               processing_level = processing_level):
                    print(f'--SIMULATION {ssp} - TROUVÉE---')
                    ds_sim_fut =pcat.search(id = id_fut, 
                                            source = pilote, 
                                            member = member,
                                            experiment = ssp,
                                             processing_level = processing_level).to_dask().load()

                    if region is not None:
                        ds_sim_fut = subset.subset_shape(ds_sim_fut, shape=region) 
                    ds_out.append(ds_sim_fut)
                else:
                    print(f'--SIMULATION {ssp} - DONNÉES MANQUANTES---')
                    ds_out.append(None)
    if months:
        print(f'On sélectionne seulement les mois {months}')
        for idx, ds in enumerate(ds_out):            
            if ds is not None:
                ds_out[idx] = ds.sel(month=months)
    return ds_out
   
def cos_weight_mean(ds, dims=None):
    try:
        weights = np.cos(np.deg2rad(ds.rlat))
        print('cos-weight weighting with rlat')
    except:
        try:
            weights = np.cos(np.deg2rad(ds.cf["latitude"]))            
        except:
            weights = np.cos(np.deg2rad(ds['lat']))
    weights.name = "weights"
    if dims is not None:
        ds_agg = ds.weighted(weights).mean(dims, keep_attrs=True)          
    else:
        ds_agg = ds.weighted(weights).mean(keep_attrs=True)      
        
    return ds_agg


def colors_ipcc(experiment):
    # IPCC WG1 style guide colors
    # AR5 https://www.ipcc.ch/site/assets/uploads/2019/04/IPCC-visual-style-guide.pdf
    # AR6 https://www.ipcc.ch/site/assets/uploads/2022/09/IPCC_AR6_WGI_VisualStyleGuide_2022.pdf
    if experiment == 'historical':
        return '#984ea3' 
    elif experiment == 'ssp126':
        return "#173C66"
    elif experiment == 'ssp245':
        return "#F79420"
    elif experiment == 'ssp370':
        return "#E71D25"
    elif experiment == 'ssp585':
        return "#951B1E"
    elif experiment == 'rcp45':
        return '#70A0CD'
    elif experiment == 'rcp85':
        return "#990002"
    else:
        #return '#984ea3'
        #return '0.5'
        return '#4daf4a'
    
    
def pcolor_map(bounds, cmap, extend):
    '''    
    Sortir un colormap pour utilisation avec pcolormesh
    bounds: list ou np.array avec les niveaux [ex.: np.arange(-30,31,5)]
    cmap: nom du colormap, ex.: 'viridis', 'RdBu_r'
    extend: ajouter des couleurs > bounds[-1] ou < bounds[0]: 'both', 'max','min', 'neither'
    '''          
    #replace first color with white
    #colors[0] = "white"
    if extend == 'max':
        # create list of 7(!) colors from colormap
        cmap = plt.cm.get_cmap(cmap,len(bounds)) 
        colors = list(cmap(np.arange(len(bounds))))
        cmap = mpl.colors.ListedColormap(colors[:-1], "")        
        # set over-color to last color of list 
        cmap.set_over(colors[-1])
        cmap.set_under('white')
        norm = mpl.colors.BoundaryNorm(bounds, ncolors=len(bounds)-1, clip=False)
    elif extend == 'both':
        cmap = plt.cm.get_cmap(cmap,len(bounds)+1) 
        colors = list(cmap(np.arange(len(bounds)+1)))
        cmap = mpl.colors.ListedColormap(colors[1:-1], "")        
        # set over-color to last color of list 
        cmap.set_over(colors[-1])
        cmap.set_under(colors[0])
        norm = mpl.colors.BoundaryNorm(bounds, ncolors=len(bounds)-1, clip=False)
    elif extend == 'neither':
        cmap = plt.cm.get_cmap(cmap)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    return cmap, norm