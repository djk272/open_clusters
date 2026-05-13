#!/usr/bin/env sh

echo "Finding Cluster Membership..."
python pm_membership.py
echo "DONE"

echo "Searching Catalogs for Cluster Members..."
python cluster_catalog_search.py
echo "DONE"

echo "Matching Cluster Members with correseponding objects in Catalog search results..."
python catalog_match_cluster_members.py
echo "DONE"