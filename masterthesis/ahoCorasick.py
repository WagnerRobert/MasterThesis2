from collections import deque

__author__ = 'delur'
class Node(object):
    def __init__(self, ch, transitions, results, fail  ):
        self.ch = ch
        self.transitions = transitions
        self.results = results
        self.fail = fail

def createSearchTree(keywords):
    root = Node(None, [], [], None)
    root.fail = root
    queue = deque([root])

    for keyword in keywords:
        current_node = root
        for char in keyword:
            new_node = None
            for transition in current_node.transitions:
                if transition.ch == char:
                    new_node = transition
                    break

            if new_node is None:
                new_node = Node(char, [], [], None)
                current_node.transitions.append(new_node)
                if current_node is root:
                    new_node.fail = root

            current_node = new_node
        current_node.results.append(keyword)

    while queue:
        current_node = queue.popleft()
        for node in current_node.transitions:
            queue.append(node)
            fail_state_node = current_node.fail
            while not any(x for x in fail_state_node.transitions if node.ch == x.ch) and fail_state_node is not root:
                fail_state_node = fail_state_node.fail

            node.fail = next((x for x in fail_state_node.transitions if node.ch == x.ch and x is not node), root)

    return root

def search(tree_root, text):
    current_node = tree_root
    for i in range(len(text)):
#    for c in text:
        transition = None
        while transition is None and transition is not tree_root:
            transition = next((x for x in current_node.transitions if (x.ch == text[i] )), current_node.fail)

        current_node = transition

        if len(transition.results) > 0:
            for result in transition.results:
                yield result, i- len(result) +1
            current_node = transition.fail


def allowOneError(searchWords):
    new_word_list = []
    for word in searchWords:
        for i in range(len(word)):
            error_word = list(word)
            for char in ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T','W','Y','V']:
                error_word[i] = char
                new_word_list.append("".join(error_word))
    return new_word_list



searchWords = ["LLMKM"]
searchWords = allowOneError(searchWords)
print searchWords
tree = createSearchTree(searchWords)
print "L-LLLLLLLLL"
print "012345678901234567"
gen = search(tree, "L-LLLLLLLLL")


while True:
    try:
        blubb = gen.next()
        print blubb
    except Exception,e:
        break