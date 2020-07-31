// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Provides a map with O(log n) insertion, removal and random access, i.e. access of the value
// relative to the i-th lowest key.
// Useful for a random selection in an ordered list with variable size.
// Implemented as an augmented red-black tree.
// Code to fix RB color violation from https://www.geeksforgeeks.org/red-black-tree-set-1-introduction-2/

#ifndef DCA_UTIL_TYPE_RANDOM_ACCESS_MAP_HPP
#define DCA_UTIL_TYPE_RANDOM_ACCESS_MAP_HPP

#include <cassert>
#include <initializer_list>
#include <functional>
#include <stack>
#include <stdexcept>

#include "dca/util/containers/fixed_size_allocator.hpp"
#include "dca/util/containers/random_access_map_iterator.hpp"

namespace dca {
namespace util {

// Precondition: elements of type Key have full order.
template <class Key, class Value, std::size_t chunk_size = 256>
class RandomAccessMap {
public:
  using const_iterator = RandomAccessMapIterator<Key, Value, true>;
  using iterator = RandomAccessMapIterator<Key, Value, false>;

  RandomAccessMap() = default;
  RandomAccessMap(const std::initializer_list<std::pair<Key, Value>>& list);
  RandomAccessMap(const RandomAccessMap& rhs);
  RandomAccessMap(RandomAccessMap&& rhs);
  RandomAccessMap& operator=(const RandomAccessMap& rhs);
  RandomAccessMap& operator=(RandomAccessMap&& rhs);

  ~RandomAccessMap();

  auto begin() const -> const_iterator;
  auto end() const -> const_iterator;

  auto begin() -> iterator;
  auto end() -> iterator;

  // Insert new key, value pair if key is not already present, and returns an iterator to the node
  // and true.
  // If the key is already present, update the value and returns an iterator to the node and false.
  auto insert(const Key& key, const Value& value) noexcept -> std::pair<iterator, bool>;
  auto insert(const std::pair<Key, Value>& pair) noexcept {
    insert(pair.first, pair.second);
  }

  // Remove the node relative to key.
  // Returns: true if the key is found and removed. False if no operation is performed.
  bool erase(const Key& key) noexcept;

  // Remove the node.
  // Precondition: the node is in the map.
  void erase(iterator it);

  // Returns the iterator associated with key.
  // If the key is not in the map, returns a null iterator.
  auto findByKey(const Key& key) const noexcept -> const_iterator;
  auto findByKey(const Key& key) noexcept -> iterator;

  // Returns true if a stored key compares equal to the argument.
  bool contains(const Key& key) const noexcept;
  bool count(const Key& key) const noexcept {
    return contains(key);
  }

  // Returns an iterator relative to the 'index'-th lowest key.
  // Precondition: 0 <= index < size()
  auto findByIndex(const std::size_t index) const -> const_iterator;
  auto findByIndex(const std::size_t index) -> iterator;

  // Number of keys stored in the map.
  std::size_t size() const noexcept {
    return root_ ? root_->subtree_size : 0;
  }

  // Returns an array of ordered keys and value pairs.
  std::vector<std::pair<Key, Value>> linearize() const noexcept;

  // For testing purposes.
  bool checkConsistency() const noexcept;

private:
  using Node = details::Node<Key, Value>;
  constexpr static auto BLACK = details::BLACK;
  constexpr static auto RED = details::RED;

  void eraseImpl(Node* to_delete) noexcept;

  const Key& get_key(const Node* node) const {
    return node->data.first;
  }

  void rightRotate(Node* node);

  void leftRotate(Node* node);

  void fixRedRed(Node* node);

  void fixDoubleBlack(Node* node);

  bool isLeftChild(const Node* node) const {
    return node->parent && node->parent->left == node;
  }

  bool isRightChild(const Node* node) const {
    return node->parent && node->parent->right == node;
  }

  auto getUncle(const Node* node) const -> Node*;

  auto getSibling(const Node* node) const -> Node*;

  void updateSubtreeSize(Node* node);

  // Inverts the parent-child relationship of the two arguments.
  void moveDown(Node* node, Node* new_parent);

  // Members
  Node* root_ = nullptr;
  dca::util::details::FixedSizeAllocator<Node, chunk_size> allocator_;
};

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>::RandomAccessMap(
    const std::initializer_list<std::pair<Key, Value>>& list) {
  for (const auto& [key, val] : list)
    insert(key, val);
}

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>::~RandomAccessMap() {
  std::stack<Node*> to_delete;
  if (root_)
    to_delete.push(root_);

  while (!to_delete.empty()) {
    auto node = to_delete.top();
    to_delete.pop();

    if (node->left)
      to_delete.push(node->left);
    if (node->right)
      to_delete.push(node->right);

    allocator_.destroy(node);
  }
}

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>::RandomAccessMap(const RandomAccessMap& rhs) {
  (*this) = rhs;
}

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>::RandomAccessMap(RandomAccessMap&& rhs) {
  (*this) = std::move(rhs);
}

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>& RandomAccessMap<Key, Value, chunk_size>::operator=(
    const RandomAccessMap<Key, Value, chunk_size>& rhs) {
  if (this != &rhs) {
    *this = std::move(RandomAccessMap());  // clear content.

    for (const auto& it : rhs)
      insert(it);
  }
  return *this;
}

template <class Key, class Value, std::size_t chunk_size>
RandomAccessMap<Key, Value, chunk_size>& RandomAccessMap<Key, Value, chunk_size>::operator=(
    RandomAccessMap<Key, Value, chunk_size>&& rhs) {
  std::swap(root_, rhs.root_);
  std::swap(allocator_, rhs.allocator_);
  return *this;
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::insert(const Key& key, const Value& val) noexcept
    -> std::pair<iterator, bool> {
  if (!root_) {
    root_ = allocator_.create(key, val, nullptr);
    root_->color = BLACK;
    return {iterator(root_), true};
  }

  Node* node = root_;
  bool done = false;

  while (!done) {
    if (key == get_key(node)) {  // Key is already present. Undo changes and return.
      node->data.second = val;
      iterator return_it = iterator(node);

      node = node->parent;
      while (node) {
        --node->subtree_size;
        node = node->parent;
      }

      return {return_it, false};
    }
    ++node->subtree_size;

    if (key < get_key(node)) {
      if (node->left == nullptr) {
        node->left = allocator_.create(key, val, node);
        done = true;
      }
      node = node->left;
    }
    else {
      if (node->right == nullptr) {
        node->right = allocator_.create(key, val, node);
        done = true;
      }
      node = node->right;
    }
  }

  // Check colors
  fixRedRed(node);

  //  assert(checkConsistency());
  return {iterator(node), true};
}

template <class Key, class Value, std::size_t chunk_size>
bool RandomAccessMap<Key, Value, chunk_size>::erase(const Key& key) noexcept {
  if (!root_)
    return false;
  Node* to_delete = root_;

  // Search while updating subtree count.
  bool found = false;

  while (true) {
    if (key == get_key(to_delete)) {
      found = true;
      break;
    }

    --to_delete->subtree_size;

    if (key < get_key(to_delete)) {
      if (!to_delete->left)
        break;
      to_delete = to_delete->left;
    }
    else {
      if (!to_delete->right)
        break;
      to_delete = to_delete->right;
    }
  }

  if (!found) {  // undo change
    ++to_delete->subtree_size;
    while (to_delete->parent) {
      to_delete = to_delete->parent;
      ++to_delete->subtree_size;
    }

    return false;
  }

  if (to_delete->left != nullptr && to_delete->right != nullptr) {  // to_delete has two children.
    Node* const original = to_delete;
    --to_delete->subtree_size;
    to_delete = to_delete->right;
    while (to_delete->left) {
      --to_delete->subtree_size;
      to_delete = to_delete->left;
    }

    original->data = std::move(to_delete->data);
  }

  --to_delete->subtree_size;
  eraseImpl(to_delete);

  return true;
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::erase(iterator it) {
  Node* to_delete = it.node_;

  if (to_delete->left != nullptr && to_delete->right != nullptr) {  // to_delete has two children.
    Node* const original = to_delete;
    to_delete = to_delete->right;
    while (to_delete->left) {
      to_delete = to_delete->left;
    }

    original->data = std::move(to_delete->data);
  }

  eraseImpl(to_delete);

  // Update subtree counts.
  Node* ancestor = to_delete->parent;
  while (ancestor) {
    --ancestor->subtree_size;
    ancestor = ancestor->parent;
  }
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::eraseImpl(Node* to_delete) noexcept {
  Node* replacement = to_delete->left ? to_delete->left : to_delete->right;

  auto color = [](const Node* n) { return n ? n->color : BLACK; };
  const bool both_black = color(replacement) == BLACK && to_delete->color == BLACK;

  if (both_black) {
    fixDoubleBlack(to_delete);
  }
  else {
    auto sibling = getSibling(to_delete);
    if (sibling && !replacement)
      sibling->color = RED;
    else if (replacement)
      replacement->color = BLACK;
  }

  // delete to_delete from the tree
  Node* parent = to_delete->parent;
  if (isLeftChild(to_delete)) {
    parent->left = replacement;
  }
  else if (isRightChild(to_delete)) {
    parent->right = replacement;
  }
  if (replacement)
    replacement->parent = parent;

  // Update root if necessary.
  if (to_delete == root_) {
    root_ = replacement;
  }

  allocator_.destroy(to_delete);

  //  assert(checkConsistency());
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::fixRedRed(Node* x) {
  // if x is root color it black and return
  if (x == root_) {
    x->color = BLACK;
    return;
  }

  // initialize relatives.
  Node* parent = x->parent;
  Node* grandparent = parent->parent;
  Node* uncle = getUncle(x);

  if (parent->color != BLACK) {
    if (uncle && uncle->color == RED) {
      // uncle is red, perform recoloring and recurse
      parent->color = BLACK;
      uncle->color = BLACK;
      grandparent->color = RED;
      fixRedRed(grandparent);
    }
    else {
      if (isLeftChild(parent)) {
        if (isLeftChild(x)) {
          // for left right
          std::swap(parent->color, grandparent->color);
        }
        else {
          leftRotate(parent);
          std::swap(x->color, grandparent->color);
        }
        // for left left and left right
        rightRotate(grandparent);
      }
      else {
        if (isLeftChild(x)) {
          // for right left
          rightRotate(parent);
          std::swap(x->color, grandparent->color);
        }
        else {
          std::swap(parent->color, grandparent->color);
        }

        // for right right and right left
        leftRotate(grandparent);
      }
    }
  }
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::findByIndex(const std::size_t index) -> iterator {
  if (index >= size())
    throw(std::out_of_range("Index out of range"));

  Node* node = root_;

  std::size_t on_the_left = 0;
  while (true) {
    assert(node);

    auto new_on_the_left = on_the_left;
    if (node->left)
      new_on_the_left += node->left->subtree_size;

    if (new_on_the_left == index) {
      return iterator(node);
    }
    else if (new_on_the_left > index) {  // go left
      node = node->left;
    }
    else {  // go right
      on_the_left = new_on_the_left + 1;
      node = node->right;
    }
  }
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::findByIndex(const std::size_t index) const
    -> const_iterator {
  return const_cast<RandomAccessMap&>(*this).findByIndex(index);
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::findByKey(const Key& key) noexcept -> iterator {
  Node* node = root_;
  while (node) {
    if (get_key(node) == key)
      return iterator(node);
    else if (key < get_key(node))
      node = node->left;
    else
      node = node->right;
  }

  // Key not found.
  return iterator(nullptr);
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::findByKey(const Key& key) const noexcept
    -> const_iterator {
  // Avoid code duplication with a cast to non-const (const iterator does not allow data modification).
  return const_iterator(const_cast<RandomAccessMap&>(*this).findByKey(key));
}

template <class Key, class Value, std::size_t chunk_size>
bool RandomAccessMap<Key, Value, chunk_size>::contains(const Key& key) const noexcept {
  const Node* node = root_;
  while (node) {
    if (get_key(node) == key)
      return true;
    else if (key < get_key(node))
      node = node->left;
    else
      node = node->right;
  }

  return false;
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::fixDoubleBlack(Node* x) {
  if (x == root_) {  // Reached root
    return;
  }

  Node* sibling = getSibling(x);
  Node* parent = x->parent;

  auto has_red_child = [](Node* n) {
    return (n->left != NULL && n->left->color == RED) || (n->right != NULL && n->right->color == RED);
  };

  if (sibling == NULL) {
    // No sibiling, double black pushed up
    fixDoubleBlack(parent);
  }
  else {
    if (sibling->color == RED) {
      // Sibling red
      parent->color = RED;
      sibling->color = BLACK;
      if (isLeftChild(sibling)) {
        // left case
        rightRotate(parent);
      }
      else {
        // right case
        leftRotate(parent);
      }
      fixDoubleBlack(x);
    }
    else {
      // Sibling black
      if (has_red_child(sibling)) {
        // at least 1 red children
        if (sibling->left != NULL and sibling->left->color == RED) {
          if (isLeftChild(sibling)) {
            // left left
            sibling->left->color = sibling->color;
            sibling->color = parent->color;
            rightRotate(parent);
          }
          else {
            // right left
            sibling->left->color = parent->color;
            rightRotate(sibling);
            leftRotate(parent);
          }
        }
        else {
          if (isLeftChild(sibling)) {
            // left right
            sibling->right->color = parent->color;
            leftRotate(sibling);
            rightRotate(parent);
          }
          else {
            // right right
            sibling->right->color = sibling->color;
            sibling->color = parent->color;
            leftRotate(parent);
          }
        }
        parent->color = BLACK;
      }
      else {
        // 2 black children
        sibling->color = RED;
        if (parent->color == BLACK)
          fixDoubleBlack(parent);
        else
          parent->color = BLACK;
      }
    }
  }
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::rightRotate(Node* const node) {
  // new parent will be node's left child
  Node* new_parent = node->left;

  // update root if current node is root
  if (node == root_)
    root_ = new_parent;

  moveDown(node, new_parent);

  // connect node with new parent's right element
  node->left = new_parent->right;
  // connect new parent's right element with node
  // if it is not nullptr
  if (new_parent->right != nullptr)
    new_parent->right->parent = node;

  // connect new parent with node
  new_parent->right = node;

  updateSubtreeSize(node);
  updateSubtreeSize(new_parent);
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::leftRotate(Node* node) {
  // new parent will be node's right child
  Node* new_parent = node->right;

  // update root_ if current node is root_
  if (node == root_)
    root_ = new_parent;

  moveDown(node, new_parent);

  // connect node with new parent's left element
  node->right = new_parent->left;
  // connect new parent's left element with node
  // if it is not nullptr
  if (new_parent->left != nullptr)
    new_parent->left->parent = node;

  // connect new parent with node
  new_parent->left = node;

  updateSubtreeSize(node);
  updateSubtreeSize(new_parent);
}

template <class Key, class Value, std::size_t chunk_size>
std::vector<std::pair<Key, Value>> RandomAccessMap<Key, Value, chunk_size>::linearize() const
    noexcept {
  std::vector<std::pair<Key, Value>> result;
  result.reserve(size());

  for (const auto& it : (*this))
    result.emplace_back(it);

  return result;
}

template <class Key, class Value, std::size_t chunk_size>
bool RandomAccessMap<Key, Value, chunk_size>::checkConsistency() const noexcept {
  bool child_parent_violation = false;
  bool red_red_violation = false;
  bool black_count_violation = false;
  bool subtree_size_violation = false;

  // Returns size of subtree.
  std::function<std::size_t(const Node*)> subtree_size = [&](const Node* node) -> std::size_t {
    if (!node)
      return 0;
    return 1 + subtree_size(node->left) + subtree_size(node->right);
  };

  // Check node consistency and returns number of black nodes in [node, leaves].
  std::function<int(const Node*)> check = [&](const Node* node) {
    if (node == nullptr)
      return 1;

    // Check parent-child relationship.
    if (node->left && node->left->parent != node)
      child_parent_violation = true;
    if (node->right && node->right->parent != node)
      child_parent_violation = true;

    // Check subtree size
    if (node->subtree_size != subtree_size(node))
      subtree_size_violation = true;

    // Check double red
    auto color = [&](const Node* n) { return n ? n->color : BLACK; };
    if (node->color == RED && (color(node->left) == RED || color(node->right) == RED))
      red_red_violation = true;

    // Check black count
    int count_left = check(node->left);
    int count_right = check(node->right);

    if (count_left != count_right)
      black_count_violation = true;

    return count_left + node->color == BLACK ? 1 : 0;
  };

  check(root_);

  return !black_count_violation && !red_red_violation && !child_parent_violation &&
         !subtree_size_violation;
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::getUncle(const Node* node) const -> Node* {
  const Node* parent = node->parent;
  if (isLeftChild(parent))
    return parent->parent->right;
  else if (isRightChild(parent))
    return parent->parent->left;
  else
    return nullptr;
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::getSibling(const Node* node) const -> Node* {
  if (isLeftChild(node))
    return node->parent->right;
  else if (isRightChild(node))
    return node->parent->left;
  else
    return nullptr;
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::moveDown(Node* node, Node* new_parent) {
  auto& parent = node->parent;
  if (isLeftChild(node)) {
    parent->left = new_parent;
  }
  else if (isRightChild(node)) {
    parent->right = new_parent;
  }
  new_parent->parent = parent;
  parent = new_parent;
}

template <class Key, class Value, std::size_t chunk_size>
void RandomAccessMap<Key, Value, chunk_size>::updateSubtreeSize(Node* node) {
  node->subtree_size = 1;
  if (node->left)
    node->subtree_size += node->left->subtree_size;
  if (node->right)
    node->subtree_size += node->right->subtree_size;
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::begin() -> iterator {
  Node* node = root_;
  while (node->left)
    node = node->left;

  return iterator(node);
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::begin() const -> const_iterator {
  return const_cast<RandomAccessMap&>(*this).begin();
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::end() const -> const_iterator {
  return const_iterator{nullptr};
}

template <class Key, class Value, std::size_t chunk_size>
auto RandomAccessMap<Key, Value, chunk_size>::end() -> iterator {
  return iterator{nullptr};
}

}  // namespace util
}  // namespace dca

#endif  // #define DCA_UTIL_TYPE_RANDOM_ACCESS_MAP_HPP
