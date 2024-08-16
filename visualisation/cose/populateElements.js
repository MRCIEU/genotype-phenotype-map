jsonInput = {
    gene: {
        id: 'IL6',
        display_name: "IL6"
    },
    studies: [
        {
            id: 'phenotype_3',
            data_type: 'phenotype',
            display_name: 'phenotype_3',
            full_name: 'phenotype_3',
            coloc_group_id: 'coloc_group_45',
            sub_studies: [
                {
                    unique_study_id: '123456'
                },
                {
                    unique_study_id: '654321'
                }
            ]
        },
        {
            id: 'phenotype_4',
            data_type: 'phenotype',
            display_name: 'phenotype_4',
            full_name: 'phenotype_4',
            coloc_group_id: 'coloc_group_45',
            sub_studies: [
                {
                    unique_study_id: '123456'
                },
                {
                    unique_study_id: '654321'
                }
            ]
        },
        {
            id: 'phenotype_5',
            data_type: 'phenotype',
            display_name: 'phenotype_5',
            full_name: 'phenotype_5',
            coloc_group_id: 'coloc_group_45',
            sub_studies: [
                {
                    unique_study_id: '123456'
                },
                {
                    unique_study_id: '654321'
                }
            ]
        },
        {
            id: 'phenotype_1',
            data_type: 'phenotype',
            display_name: 'phenotype_1',
            full_name: 'phenotype_1',
            coloc_group_id: 'coloc_group_23',
            sub_studies: [
                {
                    unique_study_id: '123456'
                },
                {
                    unique_study_id: '654321'
                }
            ]
        },
        {
            id: 'phenotype_2',
            data_type: 'phenotype',
            display_name: 'phenotype_2',
            full_name: 'phenotype_2',
            coloc_group_id: 'coloc_group_23',
            sub_studies: [
                {
                    unique_study_id: '123456'
                },
                {
                    unique_study_id: '654321'
                }
            ]
        },
        {
            id: 'eqtl_1',
            data_type: 'eqtl',
            display_name: 'eqtl_1',
            full_name: 'eqtl_1',
            tissue: 'Adipose Subcutaneous',
            coloc_group_id: 'coloc_group_23',
            sub_studies: [ { unique_study_id: '123456' } ]
        },
        {
            id: 'eqtl_2',
            data_type: 'eqtl',
            display_name: 'eqtl_2',
            full_name: 'eqtl_2',
            tissue: 'Liver',
            coloc_group_id: 'coloc_group_45',
            sub_studies: [ { unique_study_id: '123456' } ]
        },
        {
            id: 'sqtl_1',
            data_type: 'sqtl',
            display_name: 'sqtl_1',
            full_name: 'sqtl_1',
            tissue: 'Liver',
            coloc_group_id: 'coloc_group_45',
            sub_studies: [ { unique_study_id: '123456' } ]
        },
        {
            id: 'sqtl_2',
            data_type: 'sqtl',
            display_name: 'sqtl_2',
            full_name: 'sqtl_2',
            tissue: 'Adipose Subcutaneous',
            coloc_group_id: 'coloc_group_23',
            sub_studies: [ { unique_study_id: '123456' } ]
        },
    ],
    links: [
        {
            from: 'sqtl_2',
            from_data_type: 'sqtl',
            to: 'eqtl_1',
            posterior_prob: 0.94,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'eqtl_1',
            from_data_type: 'eqtl',
            to: 'phenotype_1',
            posterior_prob: 0.94,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'eqtl_1',
            from_data_type: 'eqtl',
            to: 'phenotype_2',
            posterior_prob: 0.94,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'phenotype_1',
            from_data_type: 'phenotype',
            to: 'phenotype_2',
            posterior_prob: 0.94,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'sqtl_2',
            from_data_type: 'sqtl',
            to: 'phenotype_2',
            posterior_prob: 0.94,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },




        {
            from: 'sqtl_1',
            from_data_type: 'sqtl',
            to: 'eqtl_2',
            posterior_prob: 0.87,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'eqtl_2',
            from_data_type: 'eqtl',
            to: 'phenotype_3',
            posterior_prob: 0.87,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'eqtl_2',
            from_data_type: 'eqtl',
            to: 'phenotype_4',
            posterior_prob: 0.87,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'eqtl_2',
            from_data_type: 'eqtl',
            to: 'phenotype_5',
            posterior_prob: 0.87,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'sqtl_1',
            from_data_type: 'sqtl',
            to: 'phenotype_3',
            posterior_prob: 0.87,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'sqtl_1',
            from_data_type: 'sqtl',
            to: 'phenotype_4',
            posterior_prob: 0.87,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'sqtl_1',
            from_data_type: 'sqtl',
            to: 'phenotype_5',
            posterior_prob: 0.87,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'phenotype_3',
            from_data_type: 'phenotype',
            to: 'phenotype_5',
            posterior_prob: 0.87,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
        {
            from: 'phenotype_4',
            from_data_type: 'phenotype',
            to: 'phenotype_5',
            posterior_prob: 0.87,
            candidate_snp: 'rs123456',
            posterior_explained_by_snp: 0.67
        },
    ]
}

let cy = window.cy = cytoscape({
    container: document.getElementById('cy'),
    boxSelectionEnabled: false,
    autounselectify: true,
    idealEdgeLength: 1,
    // animate: true,
    // animationDuration: 500,

    layout: {
        name: 'cose',
    },
    style: [
        {
            selector: ':parent',
            style: {
                'background-opacity': 0.05,
                'border-width': 0.3,
                'display': 'none'
            }
        },
        {
            selector: '.gene',
            style: {
                'shape': 'round-rectangle',
                'height': 20,
                'width': 40,
                "text-valign": "center",
                "text-halign": "center",
                'background-color': 'data(background)',
                'label': "data(label)",
                'font-size': '8px',
            }
        },
        {
            selector: '.tissue',
            style: {
                'background-color': '#149543',
                'shape': 'round-rectangle',
                'height': 10,
                'width': 40,
                "text-valign": "center",
                "text-halign": "center",
                'label': "data(label)",
                'font-size': '8px'
            }
        },
        {
            selector: '.result',
            style: {
                'shape': 'round-rectangle',
                'height': 10,
                'width': 40,
                "text-valign": "center",
                "text-halign": "center",
                'background-color': 'data(background)',
                'label': "data(label)",
                'font-size': '8px',
                'display': 'none'
            },
        },
        {
            selector: 'edge',
            style: {
                'curve-style': 'haystack',
                'haystack-radius': 0,
                'width': 0.7,
                'opacity': 0.3,
                'line-color': '#083a45'
            }
        },
        {
            selector: '.hidden-edge',
            style: {
                'opacity': 0,
                'line-color': '#df3a45'
            }
        }
    ],
});

nodeStyleMap = {
    gene: {
        background: '#51a2f5'
    },
    eqtl: {
        background: '#9df9ef'
    },
    mqtl: {
        background: '#edf7e6'
    },
    sqtl: {
        background: '#ffa8B6'
    },
    phenotype: {
        background: '#c2a0a9'
    }
}

//creating top level genenode
cy.add({
    data: {
        id: jsonInput.gene.id,
        level: 2,
        data_type: 'gene',
        label: jsonInput.gene.display_name,
        background: nodeStyleMap.gene.background
    },
    classes: 'gene'
})

//creating coloc groups as parent nodes, so we know how to cluster them
const unique_coloc_groups = [...new Set(jsonInput.studies.map(item => item.coloc_group_id))].filter(Boolean);
unique_coloc_groups.forEach(coloc_group => {
    cy.add({
        data: { id: coloc_group, data_type: 'coloc_group'},
        group: 'nodes'
    })
})

//creating tissues as nodes, maybe in a way to switch between views?
const unique_tissues= [...new Set(jsonInput.studies.map(item => item.tissue))].filter(Boolean);
unique_tissues.forEach(tissue => {
    cy.add({
        data: {
            id: tissue,
            data_type: 'tissue',
            level: 1,
            label: tissue,
        },
        classes: 'tissue',
        group: 'nodes'
    })
    cy.add({
        data: { source: jsonInput.gene.id, target: tissue},
        group: 'edges'
    })
})

//creating all nodes for representation
jsonInput.studies.forEach(study => {
    cy.add({
        data: {
            id: study.id,
            data_type: study.data_type,
            tissue: study.tissue,
            label: study.display_name,
            parent: study.coloc_group_id,
            sub_studies: study.sub_studies,
            background: nodeStyleMap[study.data_type].background
        },
        classes: 'result',
        group: 'nodes'
    })
})

//creating edges in graph
jsonInput.links.forEach(link => {
    if (link.from_data_type !== 'phenotype') {
        cy.add({
            data: { source: jsonInput.gene.id, target: link.from },
            group: 'edges'
        })
        //adds reverse link to and from qtl sources, so we can find all successors more easily
        cy.add({
            data: { source: link.to, target: link.from },
            group: 'edges'
        })
    }
    cy.add({
        data: { source: link.from, target: link.to },
        classes: 'hidden-edge',
        group: 'edges'
    })
})

cy.layout({ name: 'concentric' }).run();

cy.on('tap', 'node', function(event) {
    let data = event.target.data()
    if (data.data_type === 'gene' ) {
        cy.elements('node[data_type != "tissue"]').style('display', 'none')
        cy.elements('node[data_type = "tissue"]').style('display', 'element')
        cy.elements('node[data_type = "gene"]').style('display', 'element')
        cy.layout({ name: 'concentric' }).run();
    }
    else if (data.data_type === 'tissue') {
        cy.batch(function() {
            let select = 'node[tissue = "' + data.label + '"]'
            cy.elements(select).successors().style('display', 'element')
            cy.elements(select).parent().style('display', 'element')
            // cy.elements('node[data_type != "tissue"]').style('display', 'element')
            cy.elements('node[data_type = "tissue"]').style('display', 'none')
        })
        cy.layout({ name: 'cose' }).run();
    }
})
