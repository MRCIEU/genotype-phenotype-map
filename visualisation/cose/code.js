let elements = [
    {
        data: {
            id: "GENE",
            level: 3,
            label: "GENE",
            background: '#f0090c'
        },
        group: 'nodes'
    },
    {
        data: {
            id: "parent_group_1",
        },
        group: 'nodes'
    },
    {
        data: {
            id: "parent_group_2",
        },
        group: 'nodes'
    },
    {
        data: {
            id: "sQTL_1",
            label: "sQTL_1",
            parent: 'parent_group_2',
            background: '#00c9bc'
        },
        group: 'nodes'
    },
    {
        data: {
            id: "sQTL_2",
            label: "sQTL_2",
            parent: 'parent_group_1',
            background: '#00c9bc'
        },
        group: 'nodes'
    },
    {
        data: {
            id: "eQTL_1",
            label: "eQTL_1",
            parent: 'parent_group_1',
            level: 2,
            background: '#00c9bc'
        },
        group: 'nodes'
    },
    {
        data: {
            id: "eQTL_2",
            label: "eQTL_2",
            parent: 'parent_group_2',
            level: 2,
            background: '#00c9bc'
        },
        group: 'nodes'
    },
    {
        data: {
            id: "phenotype_1",
            label: "phenotype_1",
            parent: 'parent_group_1',
            background: '#ffc9bc'
        },
        group: 'nodes'
    },
    {
        data: {
            id: "phenotype_2",
            label: "phenotype_2",
            parent: 'parent_group_1',
            background: '#ffc9bc'
        },
        group: 'nodes'
    },
    {
        data: {
            id: "phenotype_3",
            label: "phenotype_3",
            parent: 'parent_group_2',
            level: 1,
            background: '#ffc9bc'
        },
        group: 'nodes'
    },
    {
        data: {
            id: "phenotype_4",
            label: "phenotype_4",
            parent: 'parent_group_2',
            level: 1,
            background: '#ffc9bc'
        },
        group: 'nodes'
    },
    {
        data: {
            id: "phenotype_5",
            label: "phenotype_5",
            parent: 'parent_group_2',
            level: 1,
            background: '#ffc9bc'
        },
        group: 'nodes'
    },
    {
        data: { source: 'GENE', target: 'sQTL_1' },
        group: 'edges',
    },
    {
        data: { source: 'GENE', target: 'sQTL_2' },
        group: 'edges',
    },
    {
        data: { source: 'GENE', target: 'eQTL_1' },
        group: 'edges',
    },
    {
        data: { source: 'GENE', target: 'eQTL_2' },
        group: 'edges',
    },
    {
        data: { source: 'sQTL_1', target: 'eQTL_2' },
        group: 'edges',
    },
    {
        data: { source: 'eQTL_1', target: 'sQTL_2' },
        group: 'edges',
    },
    {
        data: { source: 'eQTL_1', target: 'phenotype_1' },
        group: 'edges',
    },
    {
        data: { source: 'sQTL_2', target: 'phenotype_1' },
        group: 'edges',
    },
    {
        data: { source: 'eQTL_1', target: 'phenotype_2' },
        group: 'edges',
    },
    {
        data: { source: 'sQTL_2', target: 'phenotype_2' },
        group: 'edges',
    },
    {
        data: { source: 'phenotype_1', target: 'phenotype_2' },
        group: 'edges',
    },
    {
        data: { source: 'sQTL_1', target: 'phenotype_3' },
        group: 'edges',
    },
    {
        data: { source: 'eQTL_2', target: 'phenotype_3' },
        group: 'edges',
    },
    {
        data: { source: 'sQTL_1', target: 'phenotype_4' },
        group: 'edges',
    },
    {
        data: { source: 'phenotype_3', target: 'phenotype_4' },
        group: 'edges',
    },
    {
        data: { source: 'phenotype_3', target: 'phenotype_5' },
        group: 'edges',
    },
    {
        data: { source: 'phenotype_4', target: 'phenotype_5' },
        group: 'edges',
    },
    {
        data: { source: 'eQTL_2', target: 'phenotype_4' },
        group: 'edges',
    },
    {
        data: { source: 'sQTL_1', target: 'phenotype_5' },
        group: 'edges',
    },
    {
        data: { source: 'eQTL_2', target: 'phenotype_5' },
        group: 'edges',
    },
]

let cy = window.cy = cytoscape({
    container: document.getElementById('cy'),
    boxSelectionEnabled: false,
    autounselectify: true,
    animate: true,

    layout: {
        name: 'cose',
    },
    style: [
        {
            selector: ':parent',
            style: {
                'background-opacity': 0.333,
            }
        },
        {
            selector: 'node',
            style: {
                'height': 10,
                'width': 10,
                "text-valign": "center",
                "text-halign": "center",
                'background-color': 'data(background)',
                label: "data(label)",
                'font-size': '8px'
            },
        },
        { selector: 'edge',
            style: {
                'curve-style': 'haystack',
                'haystack-radius': 0,
                'width': 0.7,
                'opacity': 0.3,
                'line-color': '#083a45'
            }
        }
    ],
    elements: elements
});

