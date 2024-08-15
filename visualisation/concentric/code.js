let elements = [
    {
        data: {
            id: "GENE",
            level: 3,
            label: "GENE",
            background: '#f0090c'
        },
    },
    {
        data: {
            id: "sQTL_1",
            label: "sQTL_1",
            level: 2,
            background: '#00c9bc'
        }
    },
    {
        data: {
            id: "sQTL_2",
            label: "sQTL_2",
            level: 2,
            background: '#00c9bc'
        }
    },
    {
        data: {
            id: "eQTL_1",
            label: "eQTL_1",
            level: 2,
            background: '#00c9bc'
        }
    },
    {
        data: {
            id: "eQTL_2",
            label: "eQTL_2",
            level: 2,
            background: '#00c9bc'
        }
    },
    {
        data: {
            id: "phenotype_1",
            label: "phenotype_1",
            level: 1,
            background: '#ffc9bc'
        }
    },
    {
        data: {
            id: "phenotype_2",
            label: "phenotype_2",
            level: 1,
            background: '#ffc9bc'
        }
    },
    {
        data: {
            id: "phenotype_3",
            label: "phenotype_3",
            level: 1,
            background: '#ffc9bc'
        }
    },
    {
        data: {
            id: "phenotype_4",
            label: "phenotype_4",
            level: 1,
            background: '#ffc9bc'
        }
    },
    {
        data: {
            id: "phenotype_5",
            label: "phenotype_5",
            level: 1,
            background: '#ffc9bc'
        }
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
      name: 'concentric',
      concentric: function( node ){
        return node.data("level");
    },
    // concentric: function( node ){
    //   return node.degree();
    // },
    // levelWidth: function( nodes ){
    //   return 2;
    // }
  },

  style: [
    {
      selector: 'node',
      style: {
        'height': 3,
        'width': 3,
          "text-valign": "center",
          "text-halign": "center",
          'background-color': 'data(background)',
          label: "data(label)",
          'font-size': '1px'
      },
    },
      {
    "selector": ".center-center",
    "style": {
      "text-valign": "center",
      "text-halign": "center"
    }
  },

    {
      selector: 'edge',
      style: {
        'curve-style': 'bezier',
        'haystack-radius': 0,
        'width': 0.2,
        'opacity': 0.5,
        'line-color': '#a8eae5'
      }
    }
  ],

  elements: elements
});

